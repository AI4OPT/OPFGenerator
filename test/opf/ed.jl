function test_opf_pm(::Type{OPFGenerator.EconomicDispatch}, data::Dict)
    OPF = OPFGenerator.EconomicDispatch

    data["basic_network"] || error("Input data must be in basic format to test")
    G = length(data["gen"])
    E = length(data["branch"])

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
    # @test res["ptdf_iterations"] == res_pm["iterations"]

    # check that iterative ptdf produces same solution
    opf2 = OPFGenerator.build_opf(OPF, data, solver, iterative_ptdf=false)
    OPFGenerator.solve!(opf2)
    res2 = OPFGenerator.extract_result(opf2)
    @test isapprox(res["objective"], res2["objective"], atol=1e-6, rtol=1e-6)

    h5 = OPFGenerator.json2h5(OPF, res)
    @test haskey(h5, "meta")
    @test haskey(h5, "primal")
    @test haskey(h5, "dual")
    Gs = [
        h5["primal"]["pg"], h5["primal"]["r"],
        h5["dual"]["mu_pg_lb"], h5["dual"]["mu_pg_ub"],
        h5["dual"]["mu_r_lb"], h5["dual"]["mu_r_ub"],
        h5["dual"]["mu_total_generation"],
         # TODO: move reserve bounds to input
        h5["primal"]["rmin"], h5["primal"]["rmax"],
    ]
    Es = [
        h5["primal"]["pf"], h5["primal"]["df"],
        h5["dual"]["mu_pf_lb"], h5["dual"]["mu_pf_ub"],
        h5["dual"]["mu_df_lb"],
        h5["dual"]["lam_ptdf"],
    ]
    singles = [
        h5["primal"]["dpb_surplus"],
        h5["primal"]["dpb_shortage"],
        h5["primal"]["dr_shortage"],
        h5["dual"]["mu_dpb_surplus_lb"],
        h5["dual"]["mu_dpb_shortage_lb"],
        h5["dual"]["mu_dr_shortage_lb"],
        h5["dual"]["mu_power_balance"],
        h5["dual"]["mu_reserve_requirement"],
    ]
    @test all(size(v) == (G,) for v in Gs)
    @test all(size(v) == (E,) for v in Es)
    @test all(size(v) == () for v in singles)

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