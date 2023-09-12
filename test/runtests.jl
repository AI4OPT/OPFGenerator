using Test

using Clarabel
using Ipopt
using JuMP

using PowerModels
const PM = PowerModels
PM.silence()
using PGLib

using OPFGenerator

const IPOPT_SOLVER = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "linear_solver" => "mumps", "print_level" => 0, "tol" => 1e-6)
const CLRBL_SOLVER = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

const OPT_SOLVERS = Dict(
    PM.ACPPowerModel        => IPOPT_SOLVER,
    PM.SOCWRPowerModel      => IPOPT_SOLVER,
    PM.SOCWRConicPowerModel => CLRBL_SOLVER,
    PM.DCPPowerModel        => CLRBL_SOLVER,
)

include("opf/opf.jl")
include("standard_form.jl")
