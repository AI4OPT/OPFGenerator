function test_opf_pm(::Type{OPFGenerator.ACOPF}, data::Dict)
    OPF = OPFGenerator.ACOPF
    # Sanity checks
    data["basic_network"] || error("Input data must be in basic format to test")
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])

    # Solve OPF with PowerModels
    solver = OPT_SOLVERS[OPF]
    res_pm = PM.solve_opf(data, PM.ACPPowerModel, solver)

    # Build and solve OPF with OPFGenerator
    solver = OPT_SOLVERS[OPF]
    opf = OPFGenerator.build_opf(OPF, data, solver)
    OPFGenerator.solve!(opf)
    res = OPFGenerator.extract_result(opf)
    @test res["opf_model"] == string(OPF)
    
    # Check that problem was solved
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
        :pg => Float64[
            get(get(sol_pm["gen"], "$g", Dict()), "pg", 0) for g in 1:G
        ],
        :qg => Float64[
            get(get(sol_pm["gen"], "$g", Dict()), "qg", 0) for g in 1:G
        ],
        :va => Float64[sol_pm["bus"]["$i"]["va"] for i in 1:N],
        :vm => Float64[sol_pm["bus"]["$i"]["vm"] for i in 1:N],
    )
    model = opf.model
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
