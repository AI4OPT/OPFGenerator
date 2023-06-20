using Random
using LinearAlgebra
using StableRNGs
using TOML

BLAS.set_num_threads(1)

using PowerModels
PowerModels.silence()
using PGLib

using OPFGenerator
using JuMP
using Ipopt
using HSL

using MathOptSymbolicAD

const LIB_COINHSL = HSL.libcoinhsl


config_file = ARGS[1]
config = TOML.parsefile(config_file)
casename = config["ref"]
name = config["name"]
export_dir = config["export_dir"]

data = make_basic_network(pglib(casename))

d = Dict{String,Any}()
d["data"] = deepcopy(data)
solver = optimizer_with_attributes(Ipopt.Optimizer,
    "hsllib" => LIB_COINHSL,
    "tol" => get(config["solver"], "tol", 1e-6),
    "max_wall_time" => get(config["solver"], "max_wall_time", 3600.0),
    "linear_solver" => get(config["solver"], "linear_solver", "mumps")
)

OPFGenerator.save_json("$(export_dir)/$(name).ref.json", d)