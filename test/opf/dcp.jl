function _test_opf_detailed(opf::OPFGenerator.OPFModel{OPF}, res::Dict, res_pm::Dict) where {OPF <: PM.DCPPowerModel}
    data = opf.data
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])
    model = opf.model
    sol_pm = res_pm["solution"]

    # Force PM solution into our model, and check that the solution is feasible
    # TODO: use JuMP.primal_feasibility_report instead
    #    (would require extracting a variable => value Dict)
    var2val_pm = Dict(
        :pg => Float64[sol_pm["gen"]["$g"]["pg"] for g in 1:G],
        :va => Float64[sol_pm["bus"]["$i"]["va"] for i in 1:N],
        :pf => Float64[sol_pm["branch"]["$e"]["pf"] for e in 1:E],
    )
    # reorder pf according to branch order in our model
    dcopf_branch_order = [key[1][1] for key in keys(model[:pf])][1:E]
    var2val_pm[:pf] = var2val_pm[:pf][dcopf_branch_order]
    var2val_pm[:pf] = vcat(var2val_pm[:pf], -var2val_pm[:pf])

    @constraint(model, var2val_pm[:pg] .<= model[:pg] .<= var2val_pm[:pg])
    @constraint(model, var2val_pm[:va] .<= model[:va] .<= var2val_pm[:va])
    @constraint(model, var2val_pm[:pf] .<= model[:pf] .<= var2val_pm[:pf])

    optimize!(model)
    @test termination_status(model) ∈ [OPTIMAL, LOCALLY_SOLVED,ALMOST_LOCALLY_SOLVED]
    @test primal_status(model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]

    return nothing
end
