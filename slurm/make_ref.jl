using TOML
using PowerModels
PowerModels.silence()
using PGLib

using OPFGenerator

config_file = ARGS[1]
config = TOML.parsefile(config_file)
case_file, case_name = OPFGenerator._get_case_info(config)
export_dir = pop!(config, "export_dir")
pop!(config, "slurm")

data = OPFGenerator.OPFData(make_basic_network(PowerModels.parse_file(case_file)))

include("../exp/sampler.jl")

d = main(data, config)
d["data"] = OPFGenerator.to_dict(data)
d["config"] = config

mkpath(export_dir)
OPFGenerator.save_json("$(export_dir)/case.json", d)