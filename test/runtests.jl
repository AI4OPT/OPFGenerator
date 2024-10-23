using Random
using Test

using HDF5
using Quadmath
using StableRNGs

using Clarabel
using Ipopt
using HiGHS
using JuMP

using PowerModels
const PM = PowerModels
PM.silence()
using PGLib

using OPFGenerator

const IPOPT_SOLVER = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "linear_solver" => "mumps", "print_level" => 1, "tol" => 1e-6)
const CLRBL_SOLVER = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "verbose" => true)
const CLRBL128_SOLVER = JuMP.optimizer_with_attributes(Clarabel.Optimizer{Float128}, "verbose" => true)
const HIGHS_SOLVER = JuMP.optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => true)

const OPT_SOLVERS = Dict(
    OPFGenerator.ACOPF             => IPOPT_SOLVER,
    OPFGenerator.SOCOPFQuad        => IPOPT_SOLVER,
    OPFGenerator.SOCOPF            => CLRBL_SOLVER,
    OPFGenerator.DCOPF             => HIGHS_SOLVER,
    OPFGenerator.EconomicDispatch  => HIGHS_SOLVER,
)


const PGLIB_CASES = ["14_ieee", "30_ieee", "57_ieee", "89_pegase", "118_ieee"]
@testset "OPFGenerator" begin
    include("io.jl")
    include("bridges.jl")
    include("sampler.jl")
    include("opf/opf.jl")
    include("postprocess.jl")
end