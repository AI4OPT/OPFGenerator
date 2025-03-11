using Random
using LinearAlgebra
using TOML
using Quadmath  # for Float128 arithmetic

BLAS.set_num_threads(1)

using PowerModels
PowerModels.silence()
using PGLib

using JuMP
using Clarabel
using Ipopt
using Mosek, MosekTools
using HiGHS

using HSL_jll
const LIB_COINHSL = HSL_jll.libhsl_path

using MathOptSymbolicAD

using PGLearn

const DEFAULT_FLOAT_PRECISION = Float32

const NAME2OPTIMIZER = Dict(
    "Clarabel128" => Clarabel.Optimizer{Float128},
    "Clarabel" => Clarabel.Optimizer{Float64},
    "Ipopt" => Ipopt.Optimizer,
    "Mosek" => Mosek.Optimizer,
    "HiGHS" => HiGHS.Optimizer,
)

# Helper function to use correct arithmetic
# The default `Float64` is over-ridden only for Clarabel
_optimizer_value_type(::Any) = Float64
_optimizer_value_type(::Type{Clarabel.Optimizer{T}}) where{T} = T
_optimizer_value_type(m::MOI.OptimizerWithAttributes) = _optimizer_value_type(m.optimizer_constructor)
_optimizer_value_type(m::JuMP.AbstractModel) = JuMP.value_type(m)

function build_and_solve_model(data, config, dataset_name)
    opf_config = config["OPF"][dataset_name]
    OPF = PGLearn.OPF2TYPE[opf_config["type"]]
    solver_config = get(opf_config, "solver", Dict())

    if solver_config["name"] == "Ipopt"
        # Make sure we provide an HSL path
        get!(solver_config, "attributes", Dict())
        get!(solver_config["attributes"], "hsllib", HSL_jll.libhsl_path)
    end

    solver = optimizer_with_attributes(NAME2OPTIMIZER[solver_config["name"]],
        get(solver_config, "attributes", Dict())...
    )
    build_kwargs = Dict(Symbol(k) => v for (k, v) in get(opf_config, "kwargs", Dict()))

    tbuild = @elapsed opf = PGLearn.build_opf(OPF, data, solver;
        T=_optimizer_value_type(solver.optimizer_constructor),
        build_kwargs...
    )

    set_silent(opf.model)
    
    PGLearn.solve!(opf)

    time_extract = @elapsed res = PGLearn.extract_result(opf)
    
    res["meta"]["build_time"] = tbuild
    res["meta"]["extract_time"] = time_extract
    
    return res
end

function main(data, config)
    d = Dict{String,Any}()

    for dataset_name in keys(config["OPF"])
        d[dataset_name] = build_and_solve_model(data, config, dataset_name)
    end

    return d
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Load config file
    config = TOML.parsefile(ARGS[1])
    # Parse seed range from CL arguments
    smin = parse(Int, ARGS[2])
    smax = parse(Int, ARGS[3])

    # Grab precision used for exporting data
    fp_string = get(config, "floating_point_type", "$(DEFAULT_FLOAT_PRECISION)")
    float_type = if lowercase(fp_string) == "float32"
        Float32
    elseif lowercase(fp_string) == "float64"
        Float64
    else
        error("Invalid floating-point type: $(fp_string); only `Float32` and `Float64` are supported.")
    end
    @info "Floating-point data will be exported in `$(float_type)`"

    # Dummy run (for pre-compilation)
    data0 = PGLearn.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))
    opf_sampler0 = PGLearn.SimpleOPFSampler(data0, config["sampler"])
    rand!(MersenneTwister(1), opf_sampler0, data0)
    main(data0, config)

    # Load reference data and setup OPF sampler
    case_file, case_name = PGLearn._get_case_info(config)
    isfile(case_file) || error("Reference case file not found: $(case_file)")
    data = PGLearn.OPFData(make_basic_network(PowerModels.parse_file(case_file)))
    opf_sampler = PGLearn.SimpleOPFSampler(data, config["sampler"])

    # Data info
    N, E, L, G = data.N, data.E, data.L, data.G

    OPFs = sort(collect(keys(config["OPF"])))
    
    # Place-holder for results. 
    # For each OPF configutation, we keep a Vector of individual h5 outputs
    # These are concatenated at the end to create one H5 file per OPF configuration
    D = Dict{String,Any}()
    D["input"] = Dict{String,Any}(
        "seed" => Int[],
        # Demand data
        "pd" => Vector{Float64}[],
        "qd" => Vector{Float64}[],
        # Global reserve requirement
        "reserve_requirement" => Float64[],
        # Generator reserves min/max levels
        "rmin" => Vector{Float64}[],
        "rmax" => Vector{Float64}[],
        # Component status
        "branch_status" => Vector{Bool}[],
        "gen_status" => Vector{Bool}[],
    )
    for dataset_name in OPFs
        D[dataset_name] = Dict{String,Any}(
            "meta" => Dict{String,Any}(),
            "primal" => Dict{String,Any}(),
            "dual" => Dict{String,Any}(),
        )
    end

    # Data generation
    @info "Generating instances for case $(case_name)\nSeed range: [$smin, $smax]\nDatasets: $OPFs"
    for s in smin:smax
        rng = MersenneTwister(s)
        tgen = @elapsed data_ = rand(rng, opf_sampler)
        tsolve = @elapsed res = main(data_, config)

        # Update input data
        push!(D["input"]["seed"], s)
        push!(D["input"]["pd"], data_.pd)
        push!(D["input"]["qd"], data_.qd)
        push!(D["input"]["rmin"], data_.rmin)
        push!(D["input"]["rmax"], data_.rmax)
        push!(D["input"]["reserve_requirement"], data_.reserve_requirement)
        push!(D["input"]["branch_status"], data_.branch_status)
        push!(D["input"]["gen_status"], data_.gen_status)

        # Add output results, one for each OPF dataset
        for dataset_name in OPFs
            h = res[dataset_name]
            d = D[dataset_name]
            for (k1, v1) in h
                for (k2, v2) in v1
                    get!(d[k1], k2, Vector{typeof(v2)}())
                    push!(d[k1][k2], v2)
                end
            end
        end
        @info "Seed $s" tgen tsolve
    end
    @info "All instances completed."

    # Tensorize everything in preparation for saving to disk
    for (k1, v1) in D["input"]
        D["input"][k1] = PGLearn.tensorize(v1)
    end
    for dataset_name in OPFs
        d = D[dataset_name]
        for (k1, v1) in d, (k2, v2) in v1
            d[k1][k2] = PGLearn.tensorize(v2)
        end
        # Track random seed in dataset meta info, to simplify for post-processing
        d["meta"]["seed"] = copy(D["input"]["seed"])
    end

    # Save to disk in separate h5 files
    for (k, v) in D
        filepath = joinpath(config["export_dir"], "res_h5", "$(case_name)_$(k)_s$(smin)-s$(smax).h5")
        mkpath(dirname(filepath))
        # Convert floating-point data before exporting
        v_ = PGLearn.convert_float_data(v, float_type)
        th5write = @elapsed PGLearn.save_h5(filepath, v_)
    end

    return nothing
end
