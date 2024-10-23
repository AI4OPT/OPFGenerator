using TOML
using PowerModels
PowerModels.silence()
using PGLib

using OPFGenerator

config_file = ARGS[1]
config = TOML.parsefile(config_file)
casename = config["ref"]
export_dir = pop!(config, "export_dir")
pop!(config, "slurm")

data = OPFGenerator.OPFData(make_basic_network(pglib(casename)))

include("../exp/sampler.jl")

d = main(data, config)

mkpath(export_dir)
OPFGenerator.save_json("$(export_dir)/case.json", d)