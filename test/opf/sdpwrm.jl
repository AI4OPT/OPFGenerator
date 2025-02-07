using LinearAlgebra
using SparseArrays

function test_opf_pm(::Type{OPFGenerator.SDPOPF}, data::Dict)
    OPF = OPFGenerator.SDPOPF

    data["basic_network"] || error("Input data must be in basic format to test")
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])

    # Solve OPF with PowerModels
    solver = OPT_SOLVERS[OPF]
    res_pm = PM.solve_opf(data, PM.SDPWRMPowerModel, solver)

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
    @test isapprox(res["meta"]["primal_objective_value"], res_pm["objective"], atol=1e-3, rtol=1e-3)
    @test isapprox(res["meta"]["primal_objective_value"], res["meta"]["dual_objective_value"], rtol=1e-6)

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
        # The diagonal elements of sol_pm["WR"]. Note that the rows and columns of sol_pm["WR"] are ordered in the order of ref[:bus].
        :w => Float64[sol_pm["bus"]["$i"]["w"] for i in 1:N]
    )
    model = opf.model
    for varname in [:pg, :qg]
        x = model[varname]
        v = var2val_pm[varname]
        fix.(x, v; force=true)
    end
    fix.(diag(model[:WR]), var2val_pm[:w]; force=true)

    optimize!(model)
    @test termination_status(model) ∈ [OPTIMAL, ALMOST_OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
    # Also check that we get the same objective value as PowerModels
    @test isapprox(objective_value(opf.model), res_pm["objective"], atol=1e-3, rtol=1e-3)

    return nothing
end

"""
    _test_sdpwrm_DualFeasibility()

Test dual feasibility of SDPWRM problem.

This test is executed on the 14 bus system.
"""
function _test_sdpwrm_DualFeasibility()
    T = Float128
    data = OPFGenerator.OPFData(make_basic_network(pglib("14_ieee")))
    solver = JuMP.optimizer_with_attributes(Clarabel.Optimizer{T},
        "verbose" => true,
        "equilibrate_enable" => false,
        "tol_gap_abs"    => 1e-14,
        "tol_gap_rel"    => 1e-14,
        "tol_feas"       => 1e-14,
        "tol_infeas_rel" => 1e-14,
        "tol_ktratio"    => 1e-14,
    )
    opf = OPFGenerator.build_opf(OPFGenerator.SDPOPF, data, solver; T=T)
    # set_silent(opf.model)
    OPFGenerator.solve!(opf)
    res = OPFGenerator.extract_result(opf)

    _test_sdpwrm_DualFeasibility(data, res)

    return nothing
end

"""
    _test_sdpwrm_DualFeasibility(data, res; atol=1e-6)

Test the dual feasibility of the Semidefinite Relaxation (SDPWRM) solution.

Tests feasibility for dual constraints associated to `WR` and `WI` variables.

# Arguments
- `data::OPFData`: OPF instance data
- `res`: Result dictionary of the SOCWR optimization
- `atol=1e-6`: The absolute tolerance for feasibility checks (default is 1e-6).
"""
function _test_sdpwrm_DualFeasibility(data::OPFGenerator.OPFData, res; atol=1e-6)
    # Grab problem data
    N = data.N
    E = data.E
    br_out = data.bus_arcs_fr
    br_in = data.bus_arcs_to
    gs, bs = data.gs, data.bs
    bus_fr, bus_to = data.bus_fr, data.bus_to
    gff, gft, gtf, gtt = data.gff, data.gft, data.gtf, data.gtt
    bff, bft, btf, btt = data.bff, data.bft, data.btf, data.btt
    δθmin, δθmax = data.dvamin, data.dvamax

    # Check dual feasibility
    λp  = res["dual"]["kcl_p"]
    λq  = res["dual"]["kcl_q"]
    λpf = res["dual"]["ohm_pf"]
    λqf = res["dual"]["ohm_qf"]
    λpt = res["dual"]["ohm_pt"]
    λqt = res["dual"]["ohm_qt"]

    s = res["dual"]["s"]
    sr = res["dual"]["sr"]
    si = res["dual"]["si"]

    μθ_lb = max.(0, res["dual"]["va_diff"])
    μθ_ub = min.(0, res["dual"]["va_diff"])

    μ_w = res["dual"]["w"]
    μ_wr = res["dual"]["wr"]
    μ_wi = res["dual"]["wi"]
    
    # Reconstruct S from s, sr, si
    S = Symmetric(sparse(
        vcat([1:N;], [N+1 : 2*N;], bus_fr, bus_to, bus_fr .+ N, bus_to .+ N, bus_fr, bus_to),
        vcat([1:N;], [N+1 : 2*N;], bus_to, bus_fr, bus_to .+ N, bus_fr .+ N, bus_to .+ N, bus_fr .+ N),
        # Symmetric() uses the upper triangular part of the matrix, but there may be branches where f_bus > t_bus (entry is below the diagonal), so we repeat sr for both directions of each branch
        vcat(s, s, sr, sr, sr, sr, si, -si),
        2*N, 2*N,
        # ignore duplicate values at the same position
        (x, y) -> x
    ))

    # Check dual constraint corresponding to `WR` variables
    AR_ff_values = [
        (
            gff[e] * λpf[e] - bff[e] * λqf[e]
        )
        for e in 1:E
    ]
    AR_tt_values = [
        (
            gtt[e] * λpt[e] - btt[e] * λqt[e]
        )
        for e in 1:E
    ]
    AR_ft_values = [
        (
            gft[e] * λpf[e] 
            + gtf[e] * λpt[e]
            - bft[e] * λqf[e] 
            - btf[e] * λqt[e]
            - tan(δθmin[e]) * μθ_lb[e]
            + tan(δθmax[e]) * μθ_ub[e]
            + μ_wr[e]
        )
        for e in 1:E
    ]
    AR = Diagonal([(-gs[i] * λp[i] + bs[i] * λq[i] + μ_w[i]) for i in 1:N]) + Symmetric(sparse(
        vcat(bus_fr, bus_to, bus_fr, bus_to),
        vcat(bus_fr, bus_to, bus_to, bus_fr),
        vcat(AR_ff_values, AR_tt_values, 1/2 * AR_ft_values, 1/2 * AR_ft_values),
        N, N
    ))
    δwr = AR + S[1:N, 1:N] + S[(N+1):(2*N), (N+1):(2*N)]
    @test norm(δwr, Inf) <= atol

    # Check dual constraint corresponding to `WI` variables
    AI_values = [
        (
            bft[e] * λpf[e] 
            - btf[e] * λpt[e]
            + gft[e] * λqf[e] 
            - gtf[e] * λqt[e]
            + μθ_lb[e]
            - μθ_ub[e]
            + μ_wi[e]
        )
        for e in 1:E
    ]
    AI = sparse(
        vcat(bus_fr, bus_to), vcat(bus_to, bus_fr), vcat(1/2 * AI_values, -1/2 * AI_values), N, N
    )
    δwi = AI + S[1:N, (N+1):(2*N)] - S[(N+1):(2*N), 1:N]
    @test norm(δwi, Inf) <= atol
    return nothing
end

function _test_sdpwrm_DualSolFormat()
    data = make_basic_network(pglib("pglib_opf_case14_ieee"))
    N = length(data["bus"])
    E = length(data["branch"])

    solver = CLRBL_SOLVER
    opf = OPFGenerator.build_opf(OPFGenerator.SDPOPF, data, solver)
    set_silent(opf.model)
    OPFGenerator.solve!(opf)

    # Check shape of dual solution
    res = OPFGenerator.extract_result(opf)

    @test Set(collect(keys(res))) == Set(["meta", "primal", "dual"])
    @test size(res["dual"]["sm_fr"]) == (E, 3)
    @test size(res["dual"]["sm_to"]) == (E, 3)
    return nothing
end

function _test_sdpwrm128(data::Dict)
    opf = OPFGenerator.build_opf(OPFGenerator.SDPOPF, data, CLRBL128_SOLVER; T=Float128)

    OPFGenerator.solve!(opf)

    res = OPFGenerator.extract_result(opf)

    return nothing
end
