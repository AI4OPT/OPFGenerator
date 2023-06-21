using TOML
using PowerModels
PowerModels.silence()
using PGLib

using OPFGenerator

config_file = ARGS[1]
config = TOML.parsefile(config_file)
casename = config["ref"]
name = config["name"]
export_dir = config["export_dir"]

data = make_basic_network(pglib(casename))

OPFGenerator.save_json("$(export_dir)/$(name).ref.json", data)