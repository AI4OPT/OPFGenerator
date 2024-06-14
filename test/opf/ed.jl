function test_opf_pm(::Type{OPFGenerator.EconomicDispatch}, data::Dict)
    OPF = OPFGenerator.EconomicDispatch

    data["basic_network"] || error("Input data must be in basic format to test")
    G = length(data["gen"])

    # Build and solve OPF with OPFGenerator
    solver = OPT_SOLVERS[OPF]
    opf = OPFGenerator.build_opf(OPF, data, solver) # NOTE: this does not consider soft versions
    OPFGenerator.solve!(opf)
    res = OPFGenerator.extract_result(opf)

    # Solve OPF with PowerModels
    # must be after OPFGenerator since it modifies `data`
    solver = OPT_SOLVERS[OPF]
    res_pm = PM.solve_opf_ptdf_branch_power_cuts!(data, solver)

    # Check that the right problem was indeed solved
    @test res["opf_model"] == string(OPF)
    @test res["termination_status"] ∈ [LOCALLY_SOLVED, OPTIMAL]
    @test res["primal_status"] == FEASIBLE_POINT
    @test res["dual_status"] == FEASIBLE_POINT
    # Check objective value against PowerModels
    @test isapprox(res["objective"], res_pm["objective"], atol=1e-6, rtol=1e-6)
    @test res["ptdf_iterations"] == res_pm["iterations"]

    # Force PM solution into our model, and check that the solution is feasible
    # TODO: use JuMP.primal_feasibility_report instead
    #    (would require extracting a variable => value Dict)
    sol_pm = res_pm["solution"]
    var2val_pm = Dict(
        :pg => Float64[sol_pm["gen"]["$g"]["pg"] for g in 1:G],
        # NOTE: PowerModels does not use `pf` variables, so we only check `pg`
    )

    @constraint(opf.model, var2val_pm[:pg] .- 1e-8 .<= opf.model[:pg] .<= var2val_pm[:pg] .+ 1e-8)

    optimize!(opf.model)
    @test termination_status(opf.model) ∈ [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(opf.model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]

    return nothing
end