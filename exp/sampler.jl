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
    "Clarabel128" => Clarabel.MOIwrapper.Optimizer{Float128},
    "Clarabel" => Clarabel.MOIwrapper.Optimizer{Float64},
    "Ipopt" => Ipopt.Optimizer,
    "Mosek" => Mosek.Optimizer,
)

# Helper function to use correct arithmetic
# The default `Float64` is over-ridden only for Clarabel
value_type(::Any) = Float64
value_type(::Type{Clarabel.MOIwrapper.Optimizer{T}}) where{T} = T
value_type(m::MOI.OptimizerWithAttributes) = value_type(m.optimizer_constructor)

function main(data, config)
    d = Dict{String,Any}()
    d["meta"] = deepcopy(config)

    # Keep track of initial data file
    d["data"] = data

    # Solve all OPF formulations
    for (dataset_name, opf_config) in config["OPF"]
        OPF = OPFGenerator.OPF2TYPE[opf_config["type"]]
        solver_config = get(opf_config, "solver", Dict())

        if solver_config["name"] == "Ipopt"
            # Make sure we provide an HSL path
            # The code below does not modify anything if the user provides an HSL path
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

        # Solve OPF model
        set_silent(opf.model)
        OPFGenerator.solve!(opf)

        tsol = @elapsed res = OPFGenerator.extract_result(opf)
        res["time_build"] = tbuild
        res["time_extract"] = tsol
        d[dataset_name] = res
    end

    # Done
    return d
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Load config file
    config = TOML.parsefile(ARGS[1])
    # Parse seed range from CL arguments
    smin = parse(Int, ARGS[2])
    smax = parse(Int, ARGS[3])

    resdir_json = joinpath(config["export_dir"], "res_json")
    mkpath(resdir_json)

    # Dummy run (for pre-compilation)
    data0 = make_basic_network(pglib("14_ieee"))
    opf_sampler0 = OPFGenerator.SimpleOPFSampler(data0, config["sampler"])
    rand(StableRNG(1), opf_sampler0)
    d0 = main(data0, config)

    # Load reference data and setup OPF sampler
    data = make_basic_network(pglib(config["ref"]))
    opf_sampler = OPFGenerator.SimpleOPFSampler(data, config["sampler"])

    OPFs = keys(config["OPF"])
    caseref = config["ref"]
    
    # Data generation
    @info "Generating instances for case $caseref\nSeed range: [$smin, $smax]\nDatasets: $OPFs"
    for s in smin:smax
        rng = StableRNG(s)
        tgen = @elapsed data_ = rand(rng, opf_sampler)
        tsolve = @elapsed d = main(data_, config)
        d["meta"]["seed"] = s
        twrite = @elapsed save_json(joinpath(resdir_json, config["ref"] * "_s$s.json.gz"), d)
        @info "Seed $s" tgen tsolve twrite
    end

    @info "All instances completed."

    return nothing
end
