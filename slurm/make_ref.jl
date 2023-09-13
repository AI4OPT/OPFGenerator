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
d = main(data, config)

OPFGenerator.save_json("$(export_dir)/$(datasetname).ref.json", d)