_test_acopf(casename::String) = _test_acopf(PM.make_basic_network(pglib(casename)))

"""
    _test_acopf(data)

Solve ACOPF problem and compare to PowerModels' solution.
"""
function _test_acopf(data::Dict)
    data["basic_network"] || error("Input data must be in basic format to test")
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])

    solver = optimizer_with_attributes(Ipopt.Optimizer,
        "print_level" => 0,
    )
    res_pm = solve_ac_opf(data, solver)
    sol_pm = res_pm["solution"]

    # build and solve ACOPF 
    acopf = ACOPFGenerator.build_acopf(data, solver)
    optimize!(acopf)
    res = ACOPFGenerator._extract_solution(acopf, data)
    sol = res["solution"]

    # Check that we get consistent results with PowerModels
    @test res["termination_status"] == LOCALLY_SOLVED
    @test isapprox(res["objective"], res_pm["objective"], atol=1e-6, rtol=1e-6)

    # Force PM solution into our model, and check that the solution is feasible
    # TODO: use JuMP.primal_feasibility_report instead
    #    (would require extracting a variable => value Dict)
    var2val_pm = Dict(
        :pg => Float64[sol_pm["gen"]["$g"]["pg"] for g in 1:G],
        :qg => Float64[sol_pm["gen"]["$g"]["qg"] for g in 1:G],
        :va => Float64[sol_pm["bus"]["$i"]["va"] for i in 1:N],
        :vm => Float64[sol_pm["bus"]["$i"]["vm"] for i in 1:N],
    )
    for varname in [:pg, :qg, :va, :vm]
        x = acopf[varname]
        v = var2val_pm[varname]
        @constraint(acopf, v .<= x .<= v)
    end

    optimize!(acopf)
    @test termination_status(acopf) ∈ [LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(acopf) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]

    return nothing
end

@testset "ACOPF" begin
    @testset "$casename" for casename in ["14_ieee", "30_ieee", "57_ieee", "89_pegase", "118_ieee"]
        _test_acopf("pglib_opf_case$(casename)")
    end
end

