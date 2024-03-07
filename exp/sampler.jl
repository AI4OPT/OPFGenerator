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

function build_models(data, config)
    caseref = config["ref"]
    OPFs = keys(config["OPF"])
    @info "Building models for $caseref\nFormulations: $OPFs"

    opf_models = Dict{String, Tuple{OPFGenerator.OPFModel{<:PowerModels.AbstractPowerModel}, Float64}}()
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

        set_silent(opf.model)
        opf_models[dataset_name] = (opf, tbuild)
        @info "Built $dataset_name in $tbuild seconds"
    end

    return opf_models
end

function main(data, opf_sampler, opf_models, smin, smax, config)
    caseref = config["ref"]
    resdir_json = joinpath(config["export_dir"], "res_json")

    @info "Generating instances for case $caseref\nSeed range: [$smin, $smax]"
    for s in smin:smax
        rng = StableRNG(s)

        tgen = @elapsed rand!(rng, opf_sampler, data)

        d = Dict{String,Any}()
        d["data"] = data
        d["meta"] = config
        d["meta"]["seed"] = s

        # Solve all OPF formulations
        ttrial = @elapsed for dataset_name in keys(opf_models)
            opf = opf_models[dataset_name][1]

            tupdate = @elapsed OPFGenerator.update!(opf, data)
            tsolve = @elapsed OPFGenerator.solve!(opf)
    
            textract = @elapsed res = OPFGenerator.extract_result(opf)
            @info "$dataset_name $(termination_status(opf.model))\nUpdated in $tupdate\nSolved in $tsolve\nExtracted in $textract"
            
            d[dataset_name] = res
            d[dataset_name]["time_build"] = opf_models[dataset_name][2]
        end

        twrite = @elapsed save_json(joinpath(resdir_json, config["ref"] * "_s$s.json.gz"), d)
        @info "Seed $s" tgen ttrial twrite
    end

    @info "All instances completed."

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Load config file
    config = TOML.parsefile(ARGS[1])
    # Parse seed range from CL arguments
    smin = parse(Int, ARGS[2])
    smax = parse(Int, ARGS[3])

    resdir_json = joinpath(config["export_dir"], "res_json")
    mkpath(resdir_json)

    # Load reference data and setup OPF sampler
    data = make_basic_network(pglib(config["ref"]))
    opf_sampler = OPFGenerator.SimpleOPFSampler(data, config["sampler"])

    opf_models = build_models(data, config)
    
    main(data, opf_sampler, opf_models, smin, smax, config)

    return nothing
end
