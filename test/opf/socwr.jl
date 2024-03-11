using LinearAlgebra

const SOCWRPowerModel = Union{PM.SOCWRPowerModel,PM.SOCWRConicPowerModel}

function _test_opf_detailed(opf::OPFModel{OPF}, res::Dict, res_pm::Dict) where {OPF <: SOCWRPowerModel}
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
        :w  => Float64[sol_pm["bus"]["$i"]["w"] for i in 1:N],
    )
    for varname in [:pg, :qg, :w]
        x = model[varname]
        v = var2val_pm[varname]
        @constraint(model, v .<= x .<= v)
    end

    optimize!(model)
    @test termination_status(model) ∈ [OPTIMAL, ALMOST_OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]

    return nothing
end

"""
    _test_socwr_DualFeasibility()

Test dual feasibility of SOCWRConic problem.

This test is executed on the 118 bus system.
"""
function _test_socwr_DualFeasibility()
    T = Float128
    data = make_basic_network(pglib("pglib_opf_case118_ieee"))
    solver = JuMP.optimizer_with_attributes(Clarabel.Optimizer{T},
        "verbose" => true,
        "equilibrate_enable" => false,
        "tol_gap_abs"    => 1e-14,
        "tol_gap_rel"    => 1e-14,
        "tol_feas"       => 1e-14,
        "tol_infeas_rel" => 1e-14,
        "tol_ktratio"    => 1e-14,
    )
    opf = OPFGenerator.build_opf(SOCWRConicPowerModel, data, solver; T=T)
    # set_silent(opf.model)
    OPFGenerator.solve!(opf)
    res = OPFGenerator.extract_result(opf)

    # Dual feasibility check
    ref = PM.build_ref(data)[:it][:pm][:nw][0]
    N = length(ref[:bus])
    E = length(ref[:branch])
    bus_loads = [
        [ref[:load][l] for l in ref[:bus_loads][i]]
        for i in 1:N
    ]
    bus_shunts = [
        [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        for i in 1:N
    ]
    gs = [sum(shunt["gs"] for shunt in bus_shunts[i]; init=0.0) for i in 1:N]
    bs = [sum(shunt["bs"] for shunt in bus_shunts[i]; init=0.0) for i in 1:N]

    # Identifying entering / exiting branches
    g = [PM.calc_branch_y(ref[:branch][e])[1] for e in 1:E]
    b = [PM.calc_branch_y(ref[:branch][e])[2] for e in 1:E]
    tr = [PM.calc_branch_t(ref[:branch][e])[1] for e in 1:E]
    ti = [PM.calc_branch_t(ref[:branch][e])[2] for e in 1:E]
    ttm  = abs2.(tr) + abs2.(ti)
    g_fr = [ref[:branch][e]["g_fr"] for e in 1:E]
    g_to = [ref[:branch][e]["g_to"] for e in 1:E]
    b_fr = [ref[:branch][e]["b_fr"] for e in 1:E]
    b_to = [ref[:branch][e]["b_to"] for e in 1:E]
    br_in  = [Tuple{Int,Int,Int}[] for _ in 1:N]  # entering branches
    br_out = [Tuple{Int,Int,Int}[] for _ in 1:N]  # existing branches
    for (e, br) in ref[:branch]
        i = br["f_bus"]
        j = br["t_bus"]
        push!(br_out[i], (e, i, j))
        push!(br_in[j], (e, i, j))
    end

    # Check dual feasibility for select buses and constraints
    λp  = -[res["solution"]["bus"]["$i"]["lam_kirchhoff_active"] for i in 1:N]
    λq  = -[res["solution"]["bus"]["$i"]["lam_kirchhoff_reactive"] for i in 1:N]
    λpf = -[res["solution"]["branch"]["$e"]["lam_ohm_active_fr"] for e in 1:E]
    λqf = -[res["solution"]["branch"]["$e"]["lam_ohm_reactive_fr"] for e in 1:E]
    λpt = -[res["solution"]["branch"]["$e"]["lam_ohm_active_to"] for e in 1:E]
    λqt = -[res["solution"]["branch"]["$e"]["lam_ohm_reactive_to"] for e in 1:E]

    ωf = [res["solution"]["branch"]["$e"]["nu_voltage_prod_soc_1"] for e in 1:E]
    ωt = [res["solution"]["branch"]["$e"]["nu_voltage_prod_soc_2"] for e in 1:E]
    ωr = [res["solution"]["branch"]["$e"]["nu_voltage_prod_soc_3"] for e in 1:E]
    ωi = [res["solution"]["branch"]["$e"]["nu_voltage_prod_soc_4"] for e in 1:E]

    μ_w = [
        res["solution"]["bus"]["$i"]["mu_w_lb"] + res["solution"]["bus"]["$i"]["mu_w_ub"]
        for i in 1:N
    ]

    # Check dual constraint corresponding to `w` variables
    for i in 1:N
        @test isapprox(
            -gs[i] * λp[i]
            + bs[i] * λq[i]
            + sum(
                (+(g[e]+g_fr[e])/ttm[e]) * λpf[e]
                + (-(b[e]+b_fr[e])/ttm[e]) * λqf[e]
                + ωf[e] / sqrt(2)
                for (e, _, _) in br_out[i];
                init=zero(Float128)
            )
            + sum(
                (g[e]+g_to[e]) * λpt[e] 
                + (-(b[e]+b_to[e])) * λqt[e] 
                + ωt[e] / sqrt(2)
                for (e, _, _) in br_in[i];
                init=zero(Float128)
            )
            + μ_w[i],
            0.0;
            atol=1e-8,
            rtol=1e-8,
        )
    end

    # Check dual constraint corresponding to `wr` variables
    δθmin = [ref[:branch][e]["angmin"] for e in 1:E]
    δθmax = [ref[:branch][e]["angmax"] for e in 1:E]
    μθ_lb = [res["solution"]["branch"]["$e"]["mu_va_diff_lb"] for e in 1:E]
    μθ_ub = [-res["solution"]["branch"]["$e"]["mu_va_diff_ub"] for e in 1:E]
    for e in 1:E
        @test isapprox(
                ((-g[e]*tr[e]+b[e]*ti[e]) / ttm[e]) * λpf[e] 
            + ((-g[e]*tr[e]-b[e]*ti[e]) / ttm[e]) * λpt[e]
            - ((-b[e]*tr[e]-g[e]*ti[e]) / ttm[e]) * λqf[e] 
            - ((-b[e]*tr[e]+g[e]*ti[e]) / ttm[e]) * λqt[e]
            - tan(δθmin[e]) * μθ_lb[e]
            + tan(δθmax[e]) * μθ_ub[e]
            + ωr[e],
            0.0;
            atol=1e-8,
            rtol=1e-8,
        )
    end

    # Check dual constraint corresponding to `wi` variables
    for e in 1:E
        @test isapprox(
              ((-b[e]*tr[e]-g[e]*ti[e]) / ttm[e]) * λpf[e] 
            - ((-b[e]*tr[e]+g[e]*ti[e]) / ttm[e]) * λpt[e]
            + ((-g[e]*tr[e]+b[e]*ti[e]) / ttm[e]) * λqf[e] 
            - ((-g[e]*tr[e]-b[e]*ti[e]) / ttm[e]) * λqt[e]
            + μθ_lb[e]
            - μθ_ub[e]
            + ωi[e],
            0.0;
            atol=1e-8,
            rtol=1e-8,
        )
    end

    return nothing
end

function _test_socwr128(data::Dict)
    opf = OPFGenerator.build_opf(PM.SOCWRConicPowerModel, data, CLRBL128_SOLVER; T=Float128)

    OPFGenerator.solve!(opf)

    res = OPFGenerator.extract_result(opf)

    return nothing
end
