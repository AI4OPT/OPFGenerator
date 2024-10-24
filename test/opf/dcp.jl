function test_opf_pm(::Type{OPFGenerator.DCOPF}, data::Dict)
    OPF = OPFGenerator.DCOPF

    data["basic_network"] || error("Input data must be in basic format to test")
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])

    # Solve OPF with PowerModels
    solver = OPT_SOLVERS[OPF]
    res_pm = PM.solve_opf(data, PM.DCPPowerModel, solver)

    # Build and solve OPF with OPFGenerator
    solver = OPT_SOLVERS[OPF]
    opf = OPFGenerator.build_opf(OPF, data, solver)
    OPFGenerator.solve!(opf)
    res = OPFGenerator.extract_result(opf)

    # Check that the right problem was indeed solved
    @test res["meta"]["formulation"] == string(OPF)
    @test res["meta"]["termination_status"] ∈ ["LOCALLY_SOLVED", "OPTIMAL"]
    @test res["meta"]["primal_status"] == "FEASIBLE_POINT"
    @test res["meta"]["dual_status"] == "FEASIBLE_POINT"
    # Check objective value against PowerModels
    @test isapprox(res["meta"]["primal_objective_value"], res_pm["objective"], atol=1e-6, rtol=1e-6)
    @test isapprox(res["meta"]["primal_objective_value"], res["meta"]["dual_objective_value"], rtol=1e-6)

    # Force PM solution into our model, and check that the solution is feasible
    # TODO: use JuMP.primal_feasibility_report instead
    #    (would require extracting a variable => value Dict)
    sol_pm = res_pm["solution"]
    var2val_pm = Dict(
        :pg => Float64[
            get(get(sol_pm["gen"], "$g", Dict()), "pg", 0) for g in 1:G
        ],
        :va => Float64[sol_pm["bus"]["$i"]["va"] for i in 1:N],
        :pf => Float64[
            get(get(sol_pm["branch"], "$e", Dict()), "pf", 0) for e in 1:E
        ],
    )
    model = opf.model

    @constraint(model, var2val_pm[:pg] .<= model[:pg] .<= var2val_pm[:pg])
    @constraint(model, var2val_pm[:va] .<= model[:va] .<= var2val_pm[:va])
    @constraint(model, var2val_pm[:pf] .<= model[:pf] .<= var2val_pm[:pf])

    optimize!(model)
    @test termination_status(model) ∈ [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]

    return nothing
end
