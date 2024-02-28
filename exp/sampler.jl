using Random
using LinearAlgebra
using StableRNGs
using TOML

BLAS.set_num_threads(1)

using PowerModels
PowerModels.silence()
using PGLib

using JuMP
using Ipopt
using Mosek, MosekTools

using HSL_jll
const LIB_COINHSL = HSL_jll.libhsl_path

using MathOptSymbolicAD

using OPFGenerator

const NAME2OPTIMIZER = Dict(
    "Ipopt" => Ipopt.Optimizer,
    "Mosek" => Mosek.Optimizer,
)


function main(data, parametric_models)
    d = Dict{String,Any}()
    # Keep track of initial data file
    d["data"] = data

    # Solve all OPF formulations
    for dataset_name in keys(parametric_models)
        opf = parametric_models[dataset_name]
        OPFGenerator.change_loads(opf, data["load"])
        OPFGenerator.solve!(opf)

        textract = @elapsed res = OPFGenerator.extract_result(opf)
        res["time_extract"] = textract
        d[dataset_name] = res
    end

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

    # Load reference data and setup OPF sampler
    data = make_basic_network(pglib(config["ref"]))
    opf_sampler = OPFGenerator.SimpleOPFSampler(data, config["sampler"])

    OPFs = keys(config["OPF"])
    caseref = config["ref"]

    @info "Building parametric models for $caseref\nDatasets: $OPFs"
    parametric_models = Dict{String, OPFGenerator.OPFModel{<:PowerModels.AbstractPowerModel}}()
    build_times = Dict{String,Float64}()
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

        tbuild = @elapsed opf = OPFGenerator.build_opf(OPF, data, solver; build_kwargs...)

        set_silent(opf.model)
        parametric_models[dataset_name] = opf
        build_times[dataset_name] = tbuild
        tsolve = @elapsed OPFGenerator.solve!(opf)
        @info "Built $dataset_name in $tbuild seconds, solved in $tsolve seconds"
    end
    
    # Data generation
    @info "Generating instances for case $caseref\nSeed range: [$smin, $smax]"
    for s in smin:smax
        rng = StableRNG(s)

        tgen = @elapsed data_ = rand(rng, opf_sampler)
        tsolve = @elapsed d = main(data_, parametric_models)

        d["meta"] = config
        d["meta"]["seed"] = s
        for dataset_name in keys(parametric_models)
            d[dataset_name]["time_build"] = build_times[dataset_name]
        end
        twrite = @elapsed save_json(joinpath(resdir_json, config["ref"] * "_s$s.json.gz"), d)
        @info "Seed $s" tgen tsolve twrite
    end

    @info "All instances completed."

    return nothing
end
