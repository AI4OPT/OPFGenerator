function test_opf_pm(::Type{PGLearn.EconomicDispatch}, data::Dict)
    OPF = PGLearn.EconomicDispatch

    data["basic_network"] || error("Input data must be in basic format to test")
    G = length(data["gen"])
    E = length(data["branch"])

    # Build and solve OPF with PGLearn
    solver = OPT_SOLVERS[OPF]
    opf = PGLearn.build_opf(OPF, data, solver) # NOTE: this does not consider soft versions
    PGLearn.solve!(opf)
    res = PGLearn.extract_result(opf)

    # Solve OPF with PowerModels
    data_pm = deepcopy(data)
    solver = OPT_SOLVERS[OPF]
    res_pm = PM.solve_opf_ptdf_branch_power_cuts!(data_pm, solver)

    # Check that the right problem was indeed solved
    @test res["meta"]["formulation"] == string(OPF)
    @test res["meta"]["termination_status"] ∈ ["LOCALLY_SOLVED", "OPTIMAL"]
    @test res["meta"]["primal_status"] == "FEASIBLE_POINT"
    @test res["meta"]["dual_status"] == "FEASIBLE_POINT"
    @test isapprox(res["meta"]["primal_objective_value"], res["meta"]["dual_objective_value"], atol=1e-6, rtol=1e-6)

    # Check objective value against PowerModels
    @test isapprox(res["meta"]["primal_objective_value"], res_pm["objective"], atol=1e-6, rtol=1e-6)
    # Force PM solution into our model, and check that the solution is feasible
    # TODO: use JuMP.primal_feasibility_report instead
    #    (would require extracting a variable => value Dict)
    sol_pm = res_pm["solution"]
    var2val_pm = Dict(
        :pg => Float64[
            get(get(sol_pm["gen"], "$g", Dict()), "pg", 0) for g in 1:G
        ],
        # NOTE: PowerModels does not use `pf` variables, so we only check `pg`
    )

    @constraint(opf.model, var2val_pm[:pg] .- 1e-8 .<= opf.model[:pg] .<= var2val_pm[:pg] .+ 1e-8)

    optimize!(opf.model)
    @test termination_status(opf.model) ∈ [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(opf.model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]


    return nothing
end
