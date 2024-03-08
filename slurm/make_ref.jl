using TOML
using PowerModels
PowerModels.silence()
using PGLib

using OPFGenerator

config_file = ARGS[1]
config = TOML.parsefile(config_file)
casename = config["ref"]
export_dir = config["export_dir"]
datasetname = splitdir(export_dir)[end]

data = make_basic_network(pglib(casename))

include("../exp/sampler.jl")

d = Dict{String,Any}()
d["data"] = data
d["meta"] = config
# d["meta"]["seed"] = "ref"

opf_models = build_models(data, config)

# Solve all OPF formulations
ttrial = @elapsed for dataset_name in keys(opf_models)
    opf = opf_models[dataset_name][1]

    OPFGenerator.solve!(opf)

    res = OPFGenerator.extract_result(opf)

    d[dataset_name] = res
    d[dataset_name]["time_build"] = opf_models[dataset_name][2]
end

OPFGenerator.save_json("$(export_dir)/$(datasetname).ref.json", d)