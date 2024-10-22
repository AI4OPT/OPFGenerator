using Random
using LinearAlgebra
using StableRNGs
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

using HSL_jll
const LIB_COINHSL = HSL_jll.libhsl_path

using MathOptSymbolicAD

using OPFGenerator

const NAME2OPTIMIZER = Dict(
    "Clarabel128" => Clarabel.Optimizer{Float128},
    "Clarabel" => Clarabel.Optimizer{Float64},
    "Ipopt" => Ipopt.Optimizer,
    "Mosek" => Mosek.Optimizer,
)

# Helper function to use correct arithmetic
# The default `Float64` is over-ridden only for Clarabel
value_type(::Any) = Float64
value_type(::Type{Clarabel.Optimizer{T}}) where{T} = T
value_type(m::MOI.OptimizerWithAttributes) = value_type(m.optimizer_constructor)

function build_and_solve_model(data, config, dataset_name)
    opf_config = config["OPF"][dataset_name]
    OPF = OPFGenerator.OPF2TYPE[opf_config["type"]]
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

    tbuild = @elapsed opf = OPFGenerator.build_opf(OPF, data, solver;
        T=value_type(solver.optimizer_constructor),
        build_kwargs...
    )

    set_silent(opf.model)
    
    OPFGenerator.solve!(opf)

    time_extract = @elapsed res = OPFGenerator.extract_result(opf)
    
    res["meta"]["build_time"] = tbuild
    res["meta"]["extract_time"] = time_extract
    
    return res
end

function main(data, config)
    d = Dict{String,Any}()
    d["meta"] = deepcopy(config)
    
    # Keep track of input data
    d["data"] = data

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

    # Dummy run (for pre-compilation)
    data0 = OPFGenerator.OPFData(make_basic_network(pglib("14_ieee")))
    opf_sampler0 = OPFGenerator.SimpleOPFSampler(data0, config["sampler"])
    rand!(StableRNG(1), opf_sampler0, data0)
    main(data0, config)

    # Load reference data and setup OPF sampler
    data = OPFGenerator.OPFData(make_basic_network(pglib(config["ref"])))
    opf_sampler = OPFGenerator.SimpleOPFSampler(data, config["sampler"])

    # Data info
    N, E, L, G = data.N, data.E, data.L, data.G

    OPFs = sort(collect(keys(config["OPF"])))
    caseref = config["ref"]
    
    # Place-holder for results. 
    # For each OPF configutation, we keep a Vector of individual h5 outputs
    # These are concatenated at the end to create one H5 file per OPF configuration
    D = Dict{String,Any}()
    D["input"] = Dict{String,Any}(
        "meta" => Dict{String,Any}("seed" => Int[]),
        "data" => Dict{String,Any}(
            "pd" => Vector{Float64}[],
            "qd" => Vector{Float64}[],
            "branch_status" => Vector{Bool}[],
            "gen_status" => Vector{Bool}[],
        )
    )
    for dataset_name in OPFs
        D[dataset_name] = Dict{String,Any}(
            "meta" => Dict{String,Any}(),
            "primal" => Dict{String,Any}(),
            "dual" => Dict{String,Any}(),
        )
    end

    # Data generation
    @info "Generating instances for case $caseref\nSeed range: [$smin, $smax]\nDatasets: $OPFs"
    for s in smin:smax
        rng = StableRNG(s)
        tgen = @elapsed data_ = rand(rng, opf_sampler)
        tsolve = @elapsed res = main(data_, config)
        res["meta"]["seed"] = s

        # Update input data
        push!(D["input"]["meta"]["seed"], s)
        push!(D["input"]["data"]["pd"], data_.pd)
        push!(D["input"]["data"]["qd"], data_.qd)
        push!(D["input"]["data"]["branch_status"], data_.branch_status)
        push!(D["input"]["data"]["gen_status"], data_.gen_status)

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
    D["input"]["meta"]["seed"] = OPFGenerator.tensorize(D["input"]["meta"]["seed"])
    for (k1, v1) in D["input"]["data"]
        D["input"]["data"][k1] = OPFGenerator.tensorize(v1)
    end
    for dataset_name in OPFs
        d = D[dataset_name]
        for (k1, v1) in d, (k2, v2) in v1
            d[k1][k2] = OPFGenerator.tensorize(v2)
        end
        # Track random seed in dataset meta info, to simplify for post-processing
        d["meta"]["seed"] = copy(D["input"]["meta"]["seed"])
    end

    # Save to disk in separate h5 files
    for (k, v) in D
        filepath = joinpath(config["export_dir"], "res_h5", "$(caseref)_$(k)_s$(smin)-s$(smax).h5")
        mkpath(dirname(filepath))
        th5write = @elapsed OPFGenerator.save_h5(filepath, v)
    end

    return nothing
end