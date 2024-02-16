using OPFGenerator: OPFModel

function test_opf(::Type{OPF}, casename::String) where {OPF <: PM.AbstractPowerModel}
    data = make_basic_network(pglib(casename))
    return test_opf(OPF, data)
end

"""
    test_opf_pm(OPF, data)

Build & solve opf using OPFGenerator, and compare against PowerModels implementation.

Returns the `OPFModel` struct, the result dictionary, and the PowerModels result dictionary.
"""
function test_opf(::Type{OPF}, data::Dict) where{OPF <: PM.AbstractPowerModel}
    # Sanity checks
    data["basic_network"] || error("Input data must be in basic format to test")
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])

    # Solve OPF with PowerModels
    solver = OPT_SOLVERS[OPF]
    res_pm = PM.solve_opf(data, OPF, solver)

    # Build and solve OPF with OPFGenerator
    opf = OPFGenerator.build_opf(OPF, data, solver)
    optimize!(opf.model)
    res = OPFGenerator.extract_result(opf)
    @test res["opf_model"] == string(OPF)
    
    # Check that problem was solved
    @test res["termination_status"] ∈ [LOCALLY_SOLVED, OPTIMAL]
    @test res["primal_status"] == FEASIBLE_POINT
    @test res["dual_status"] == FEASIBLE_POINT
    @test isapprox(res["objective"], res_pm["objective"], atol=1e-6, rtol=1e-6)

    _test_opf_detailed(opf, res, res_pm)

    return opf, res, res_pm
end

function _test_opf_detailed(opf::OPFModel{OPF}, res, res_pm) where{OPF<:PM.AbstractPowerModel}
    error("""Detailed tests not implemented for OPF formulation $(OPF)
    You must implement a function with the following signature:
        function _test_opf_detailed(opf::OPFModel{OPF}, res::Dict, res_pm::Dict) where{OPF <: $(OPF)}
            # unit tests ...
            return nothing
        end
    """)
    return nothing
end

include("acp.jl")
include("dcp.jl")
# TODO: add PTDF tests
include("socwr.jl")

# other tests
include("quad_obj.jl")

const PGLIB_CASES = ["14_ieee", "30_ieee", "57_ieee", "89_pegase", "118_ieee"]

@testset "OPF" begin
    @testset "$(OPF)" for OPF in OPFGenerator.SUPPORTED_OPF_MODELS
        @testset "$(casename)" for casename in PGLIB_CASES
            test_opf(OPF, "pglib_opf_case$(casename)")
        end

        @testset "QuadObj" begin test_quad_obj_warn(OPF) end
    end
end
