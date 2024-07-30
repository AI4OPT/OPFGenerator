using LinearAlgebra
using SparseArrays

function test_opf_pm(::Type{PM.SDPWRMPowerModel}, data::Dict)
    OPF = PM.SDPWRMPowerModel
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
    @test isapprox(res["objective"], res["objective_lb"], rtol=1e-6)

    # Force PM solution into our model, and check that the solution is feasible
    # TODO: use JuMP.primal_feasibility_report instead
    #    (would require extracting a variable => value Dict)
    sol_pm = res_pm["solution"]
    var2val_pm = Dict(
        :pg => Float64[sol_pm["gen"]["$g"]["pg"] for g in 1:G],
        :qg => Float64[sol_pm["gen"]["$g"]["qg"] for g in 1:G],
        :wm => Float64[sol_pm["bus"]["$i"]["w"] for i in 1:N]  # the diagonal elements of sol_pm["WR"]. Note that the rows and columns of sol_pm["WR"] are ordered in the order of ref[:bus].
    )
    model = opf.model
    for varname in [:pg, :qg]
        x = model[varname]
        v = var2val_pm[varname]
        @constraint(model, v .- 1e-3 .<= x .<= v .+ 1e-3)  # more tolerance for SDP
    end
    @constraint(model, var2val_pm[:wm] .- 1e-3 .<= diag(model[:WR]) .<= var2val_pm[:wm] .+ 1e-3)

    optimize!(model)
    @test termination_status(model) ∈ [OPTIMAL, ALMOST_OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
    # Also check that we get the same objective value as PowerModels
    @test isapprox(objective_value(opf.model), res_pm["objective"], atol=1e-6, rtol=1e-6)

    return nothing
end

"""
    _test_sdpwrm_DualFeasibility()

Test dual feasibility of SDPWRM problem.

This test is executed on the 118 bus system.
"""
function _test_sdpwrm_DualFeasibility()
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
    opf = OPFGenerator.build_opf(SDPWRMPowerModel, data, solver; T=T)
    # set_silent(opf.model)
    OPFGenerator.solve!(opf)
    res = OPFGenerator.extract_result(opf)

    _test_sdpwrm_DualFeasibility(data, res)

    return nothing
end

function _get_sym(i, j, v, N)
    # return an N * N sparse symmetric matrix where the (i, j) entry is v
    return Symmetric(sparse([i], [j], v, N, N))
end

function _get_skew_sym(i, j, v, N)
    # return an N * N sparse skew-symmetric matrix where the (i, j) entry is v
    return sparse([i, j], [j, i], [v, -v], N, N)
end

function _test_sdpwrm_DualFeasibility(data, res; atol=1e-6)
    ref = PM.build_ref(data)[:it][:pm][:nw][0]
    N = length(ref[:bus])
    E = length(ref[:branch])
    bus_loads = [
        [ref[:load][l] for l in ref[:bus_loads][i]]
        for i in 1:N
    ]
    # Bus-level data
    bus_shunts = [
        [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        for i in 1:N
    ]
    gs = [sum(shunt["gs"] for shunt in bus_shunts[i]; init=0.0) for i in 1:N]
    bs = [sum(shunt["bs"] for shunt in bus_shunts[i]; init=0.0) for i in 1:N]

    # Extract branch-level data
    g = [PM.calc_branch_y(ref[:branch][e])[1] for e in 1:E]
    b = [PM.calc_branch_y(ref[:branch][e])[2] for e in 1:E]
    tr = [PM.calc_branch_t(ref[:branch][e])[1] for e in 1:E]
    ti = [PM.calc_branch_t(ref[:branch][e])[2] for e in 1:E]
    ttm  = abs2.(tr) + abs2.(ti)
    g_fr = [ref[:branch][e]["g_fr"] for e in 1:E]
    g_to = [ref[:branch][e]["g_to"] for e in 1:E]
    b_fr = [ref[:branch][e]["b_fr"] for e in 1:E]
    b_to = [ref[:branch][e]["b_to"] for e in 1:E]
    δθmin = [ref[:branch][e]["angmin"] for e in 1:E]
    δθmax = [ref[:branch][e]["angmax"] for e in 1:E]

    # Check dual feasibility for select buses and constraints
    λp  = [res["solution"]["bus"]["$i"]["lam_kirchhoff_active"] for i in 1:N]
    λq  = [res["solution"]["bus"]["$i"]["lam_kirchhoff_reactive"] for i in 1:N]
    λpf = [res["solution"]["branch"]["$e"]["lam_ohm_active_fr"] for e in 1:E]
    λqf = [res["solution"]["branch"]["$e"]["lam_ohm_reactive_fr"] for e in 1:E]
    λpt = [res["solution"]["branch"]["$e"]["lam_ohm_active_to"] for e in 1:E]
    λqt = [res["solution"]["branch"]["$e"]["lam_ohm_reactive_to"] for e in 1:E]

    S = res["solution"]["S"]

    μθ_lb = [res["solution"]["branch"]["$e"]["mu_va_diff_lb"] for e in 1:E]
    μθ_ub = [-res["solution"]["branch"]["$e"]["mu_va_diff_ub"] for e in 1:E]

    μ_w = [
        res["solution"]["bus"]["$i"]["mu_w_lb"] + res["solution"]["bus"]["$i"]["mu_w_ub"]
        for i in 1:N
    ]

    # Check dual constraint corresponding to `wr` variables
    AR = (
        sum(
            (-gs[i] * λp[i] + bs[i] * λq[i]) * _get_sym(i, i, 1, N)
            for i in 1:N;
            init=0
        )
        + sum(
            λpf[e] * (
                (g[e]+g_fr[e]) / ttm[e] * _get_sym(branch["f_bus"], branch["f_bus"], 1, N)
                + (-g[e]*tr[e]+b[e]*ti[e]) / ttm[e] * _get_sym(branch["f_bus"], branch["t_bus"], 1/2, N)
            )
            + λpt[e] * (
                (g[e]+g_to[e]) * _get_sym(branch["t_bus"], branch["t_bus"], 1, N)
                + (-g[e]*tr[e]-b[e]*ti[e]) / ttm[e] * _get_sym(branch["f_bus"], branch["t_bus"], 1/2, N)
            )
            + λqf[e] * (
                -(b[e]+b_fr[e]) / ttm[e] * _get_sym(branch["f_bus"], branch["f_bus"], 1, N)
                - (-b[e]*tr[e]-g[e]*ti[e]) / ttm[e] * _get_sym(branch["f_bus"], branch["t_bus"], 1/2, N)
            )
            + λqt[e] * (
                -(b[e]+b_to[e]) * _get_sym(branch["t_bus"], branch["t_bus"], 1, N)
                - (-b[e]*tr[e]+g[e]*ti[e]) / ttm[e] * _get_sym(branch["f_bus"], branch["t_bus"], 1/2, N)
            )
            + (-tan(δθmin[e]) * μθ_lb[e] + tan(δθmax[e]) * μθ_ub[e]) * _get_sym(branch["f_bus"], branch["t_bus"], 1/2, N)
            + μ_w[i] * _get_sym(branch["f_bus"], branch["f_bus"], 1, N)
            for (e, branch) in ref[:branch];
            init=0
        )
    )
    @test norm(AR + S[1:N, 1:N] + S[(N+1):(2*N), (N+1):(2*N)], Inf) <= atol

    # Check dual constraint corresponding to `wi` variables
    AI = sum(
        _get_skew_sym(branch["f_bus"], branch["t_bus"], 1/2, N) * (
            λpf[e] * (-b[e]*tr[e]-g[e]*ti[e]) / ttm[e]
            - λpt[e] * (-b[e]*tr[e]+g[e]*ti[e]) / ttm[e]
            + λqf[e] * (-g[e]*tr[e]+b[e]*ti[e]) / ttm[e]
            - λqt[e] * (-g[e]*tr[e]-b[e]*ti[e]) / ttm[e]
            + μθ_lb[e]
            - μθ_ub[e]
        )
        for (e, branch) in ref[:branch];
        init=0
    )
    @test norm(AI + S[1:N, (N+1):(2*N)] - S[(N+1):(2*N), 1:N], Inf) <= atol
    return nothing
end

function _test_sdpwrm_DualSolFormat()
    data = make_basic_network(pglib("pglib_opf_case118_ieee"))
    N = length(data["bus"])
    E = length(data["branch"])

    solver = CLRBL_SOLVER
    opf = OPFGenerator.build_opf(SDPWRMPowerModel, data, solver)
    set_silent(opf.model)
    OPFGenerator.solve!(opf)

    # Check shape of dual solution
    res = OPFGenerator.extract_result(opf)
    @test size(res["solution"]["branch"]["1"]["nu_sm_to"]) == (3,)
    @test size(res["solution"]["branch"]["1"]["nu_sm_fr"]) == (3,)
    @test size(res["solution"]["S"]) == (2*N, 2*N)

    # Check conversion to H5 format
    h5 = OPFGenerator.json2h5(SDPWRMPowerModel, res)

    @test Set(collect(keys(h5))) == Set(["meta", "primal", "dual"])
    @test size(h5["dual"]["nu_sm_fr"]) == (E, 3)
    @test size(h5["dual"]["nu_sm_to"]) == (E, 3)
    @test size(h5["dual"]["S"]) == (2*N, 2*N)
    return nothing
end

function _test_sdpwrm128(data::Dict)
    opf = OPFGenerator.build_opf(PM.SDPWRMPowerModel, data, CLRBL128_SOLVER; T=Float128)

    OPFGenerator.solve!(opf)

    res = OPFGenerator.extract_result(opf)

    return nothing
end