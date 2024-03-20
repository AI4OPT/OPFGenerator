function test_opf_pm(::Type{PM.DCPPowerModel}, data::Dict)
    OPF = PM.DCPPowerModel

    data["basic_network"] || error("Input data must be in basic format to test")
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])

    # Solve OPF with PowerModels
    solver = OPT_SOLVERS[OPF]
    res_pm = PM.solve_opf(data, OPF, solver)

    # Build and solve OPF with OPFGenerator
    solver = OPT_SOLVERS[OPF]
    opf = OPFGenerator.build_opf(OPF, data, solver)
    OPFGenerator.solve!(opf)
    res = OPFGenerator.extract_result(opf)

    # Check that the right problem was indeed solved
    @test res["opf_model"] == string(OPF)
    @test res["termination_status"] ∈ [LOCALLY_SOLVED, OPTIMAL]
    @test res["primal_status"] == FEASIBLE_POINT
    @test res["dual_status"] == FEASIBLE_POINT
    # Check objective value against PowerModels
    @test isapprox(res["objective"], res_pm["objective"], atol=1e-6, rtol=1e-6)

    # Force PM solution into our model, and check that the solution is feasible
    # TODO: use JuMP.primal_feasibility_report instead
    #    (would require extracting a variable => value Dict)
    sol_pm = res_pm["solution"]
    var2val_pm = Dict(
        :pg => Float64[sol_pm["gen"]["$g"]["pg"] for g in 1:G],
        :va => Float64[sol_pm["bus"]["$i"]["va"] for i in 1:N],
        :pf => Float64[sol_pm["branch"]["$e"]["pf"] for e in 1:E],
    )
    model = opf.model
    # reorder pf according to branch order in our model
    dcopf_branch_order = [key[1][1] for key in keys(model[:pf])][1:E]
    var2val_pm[:pf] = var2val_pm[:pf][dcopf_branch_order]
    var2val_pm[:pf] = vcat(var2val_pm[:pf], -var2val_pm[:pf])

    @constraint(model, var2val_pm[:pg] .<= model[:pg] .<= var2val_pm[:pg])
    @constraint(model, var2val_pm[:va] .<= model[:va] .<= var2val_pm[:va])
    @constraint(model, var2val_pm[:pf] .<= model[:pf] .<= var2val_pm[:pf])

    optimize!(model)
    @test termination_status(model) ∈ [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]

    return nothing
end
