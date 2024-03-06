using Random
using Test

using Quadmath
using StableRNGs

using Clarabel
using Ipopt
using JuMP

using PowerModels
const PM = PowerModels
PM.silence()
using PGLib

using OPFGenerator

const IPOPT_SOLVER = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "linear_solver" => "mumps", "print_level" => 1, "tol" => 1e-6)
const CLRBL_SOLVER = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "verbose" => true)
const CLRBL128_SOLVER = JuMP.optimizer_with_attributes(Clarabel.Optimizer{Float128},
    "verbose" => true,
)

const OPT_SOLVERS = Dict(
    PM.ACPPowerModel               => IPOPT_SOLVER,
    PM.SOCWRPowerModel             => IPOPT_SOLVER,
    PM.SOCWRConicPowerModel        => CLRBL_SOLVER,
    PM.DCPPowerModel               => CLRBL_SOLVER,
)

include("io.jl")
include("sampler.jl")
include("opf/opf.jl")
