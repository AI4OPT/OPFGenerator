_test_soc_opf(casename::String) = _test_soc_opf(PM.make_basic_network(pglib(casename)))

"""
    _test_soc_opf(data)

Solve ACOPF problem and compare to PowerModels' solution.
"""
function _test_soc_opf(data::Dict)
    data["basic_network"] || error("Input data must be in basic format to test")
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])

    solver = optimizer_with_attributes(Ipopt.Optimizer,
        "print_level" => 0,
    )
    res_pm = solve_opf(data, SOCWRPowerModel, solver)
    sol_pm = res_pm["solution"]

    # build and solve SOC-OPF 
    socp = OPFGenerator.build_soc_opf(data, solver)
    optimize!(socp)
    res = OPFGenerator._extract_solution(socp, data)
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
        :w  => Float64[sol_pm["bus"]["$i"]["w"] for i in 1:N],
    )
    for varname in [:pg, :qg, :w]
        x = socp[varname]
        v = var2val_pm[varname]
        @constraint(socp, v .<= x .<= v)
    end

    optimize!(socp)
    @test termination_status(socp) ∈ [LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(socp) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]

    return nothing
end

@testset "SOC-OPF" begin
    @testset "$casename" for casename in ["14_ieee", "30_ieee", "57_ieee", "89_pegase", "118_ieee"]
        _test_soc_opf("pglib_opf_case$(casename)")
    end
end

