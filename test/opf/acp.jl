function _test_opf_detailed(opf::OPFGenerator.OPFModel{OPF}, res::Dict, res_pm::Dict) where {OPF <: PM.ACPPowerModel}
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
        :qg => Float64[sol_pm["gen"]["$g"]["qg"] for g in 1:G],
        :va => Float64[sol_pm["bus"]["$i"]["va"] for i in 1:N],
        :vm => Float64[sol_pm["bus"]["$i"]["vm"] for i in 1:N],
    )
    for varname in [:pg, :qg, :va, :vm]
        x = model[varname]
        v = var2val_pm[varname]
        # Ipopt will complain if we fix too many variables
        # To avoid this, we add a small tolerance below
        @constraint(model, v .- 1e-8 .<= x .<= v .+ 1e-8)
    end

    optimize!(model)
    @test termination_status(model) ∈ [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]

    return nothing
end
