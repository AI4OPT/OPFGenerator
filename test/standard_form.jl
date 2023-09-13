_test_standard_form(casename::String) = _test_standard_form(PM.make_basic_network(pglib(casename)))

"""
    _test_standard_form(data)

Convert DCOPF to standard form and compare solutions.
"""
function _test_standard_form(data::Dict)
    data["basic_network"] || error("Input data must be in basic format to test")

    # solve dcopf like usual
    solver = OPT_SOLVERS[PM.DCPPowerModel]
    dcopf = OPFGenerator.build_opf(PM.DCPPowerModel, data, solver)
    optimize!(dcopf.model)
    @test termination_status(dcopf.model) ∈ [OPTIMAL, LOCALLY_SOLVED,ALMOST_LOCALLY_SOLVED]
    @test primal_status(dcopf.model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
    dcopf_sol = Dict(ref=>value.(ref) for ref in all_variables(dcopf.model))

    # solve converted dcopf
    std_model, std = OPFGenerator.make_standard_form(dcopf, solver)
    optimize!(std_model)
    @test termination_status(std_model) ∈ [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(std_model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]

    # make sure objective values are the same
    @test isapprox(objective_value(dcopf.model), objective_value(std_model), atol=1e-6, rtol=1e-6)
    std_sol = OPFGenerator.map_standard_form_solution(std_model, std)

    # test that solutions are the same
    dcopf_vars = collect(keys(dcopf_sol))
    @test all(isapprox.(dcopf_sol[dcopf_var], std_sol[dcopf_var], atol=1e-6, rtol=1e-6) for dcopf_var in dcopf_vars)

    # force standard form solution into dcopf and see if its feasible
    for dcopf_var in dcopf_vars
        @constraint(dcopf.model, std_sol[dcopf_var] <= dcopf_var <= std_sol[dcopf_var])
    end
    optimize!(dcopf.model)
    @test termination_status(dcopf.model) ∈ [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(dcopf.model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]

    return nothing
end

@testset "StandardForm" begin
    @testset "$casename" for casename in PGLIB_CASES
        _test_standard_form("pglib_opf_case$(casename)")
    end
end
