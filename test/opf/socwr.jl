using LinearAlgebra

function test_opf_pm(::Type{OPF}, data::Dict) where {OPF <: Union{PGLearn.SOCOPFQuad,PGLearn.SOCOPF}}
    data["basic_network"] || error("Input data must be in basic format to test")
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])

    pm_type = OPF == PGLearn.SOCOPF ? PM.SOCWRConicPowerModel : PM.SOCWRPowerModel

    # Solve OPF with PowerModels
    solver = OPT_SOLVERS[OPF]
    res_pm = PM.solve_opf(data, pm_type, solver)

    # Build and solve OPF with PGLearn
    solver = OPT_SOLVERS[OPF]
    opf = PGLearn.build_opf(OPF, data, solver)
    PGLearn.solve!(opf)
    res = PGLearn.extract_result(opf)

    # Check that the right problem was indeed solved
    @test res["meta"]["formulation"] == string(OPF)
    @test res["meta"]["termination_status"] ∈ ["LOCALLY_SOLVED", "OPTIMAL"]
    @test res["meta"]["primal_status"] == "FEASIBLE_POINT"
    @test res["meta"]["dual_status"] == "FEASIBLE_POINT"
    # ⚠ we do not check against PowerModels' objective value, 
    #   because our SOC formulation is not equivalent
    # Check that primal/dual objectives are matching only for conic form
    #   (Ipopt is not good with dual objective value)
    if OPF == PGLearn.SOCOPF
        @test isapprox(res["meta"]["primal_objective_value"], res["meta"]["dual_objective_value"], rtol=1e-6)
    end

    # Force PM solution into our model, and check that the solution is feasible
    # TODO: use JuMP.primal_feasibility_report instead
    #    (would require extracting a variable => value Dict)
    # PowerModels' SOCWR formulation is more restricted than ours,
    #   so the PowerModels primal solution should be feasible
    sol_pm = res_pm["solution"]
    var2val_pm = Dict(
        :pg => Float64[
            get(get(sol_pm["gen"], "$g", Dict()), "pg", 0) for g in 1:G
        ],
        :qg => Float64[
            get(get(sol_pm["gen"], "$g", Dict()), "qg", 0) for g in 1:G
        ],
        :w  => Float64[sol_pm["bus"]["$i"]["w"] for i in 1:N],
    )
    model = opf.model
    for varname in [:pg, :qg, :w]
        x = model[varname]
        v = var2val_pm[varname]
        @constraint(model, v .<= x .<= v)
    end

    optimize!(model)
    @test termination_status(model) ∈ [OPTIMAL, ALMOST_OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
    # Also check that we get the same objective value as PowerModels
    @test isapprox(objective_value(opf.model), res_pm["objective"], atol=1e-6, rtol=1e-6)

    return nothing
end

"""
    _test_socwr_DualFeasibility()

Test dual feasibility of SOCWRConic problem.

This test is executed on the 118 bus system.
"""
function _test_socwr_DualFeasibility()
    T = Float128
    data = PGLearn.OPFData(make_basic_network(pglib("pglib_opf_case118_ieee")))
    solver = JuMP.optimizer_with_attributes(Clarabel.Optimizer{T},
        "verbose" => true,
        "equilibrate_enable" => false,
        "tol_gap_abs"    => 1e-14,
        "tol_gap_rel"    => 1e-14,
        "tol_feas"       => 1e-14,
        "tol_infeas_rel" => 1e-14,
        "tol_ktratio"    => 1e-14,
    )
    opf = PGLearn.build_opf(PGLearn.SOCOPF, data, solver; T=T)
    # set_silent(opf.model)
    PGLearn.solve!(opf)
    res = PGLearn.extract_result(opf)

    _test_socwr_DualFeasibility(data, res)

    return nothing
end


"""
    _test_socwr_DualFeasibility(data, res; atol=1e-6)

Test the dual feasibility of the Second-Order Cone Relaxation (SOCWR) solution.

Tests feasibility for dual constraints associated to `w`, `wr`, and `wi` variables.

# Arguments
- `data::OPFData`: OPF instance data
- `res`: Result dictionary of the SOCWR optimization
- `atol=1e-6`: The absolute tolerance for feasibility checks (default is 1e-6).
"""
function _test_socwr_DualFeasibility(data::PGLearn.OPFData, res; atol=1e-6)
    # Grab problem data
    N = data.N
    E = data.E
    br_out = data.bus_arcs_fr
    br_in = data.bus_arcs_to
    gs, bs = data.gs, data.bs
    gff, gft, gtf, gtt = data.gff, data.gft, data.gtf, data.gtt
    bff, bft, btf, btt = data.bff, data.bft, data.btf, data.btt
    δθmin, δθmax = data.dvamin, data.dvamax

    # Check dual feasibility for select buses and constraints
    λp  = res["dual"]["kcl_p"]
    λq  = res["dual"]["kcl_q"]
    λpf = res["dual"]["ohm_pf"]
    λqf = res["dual"]["ohm_qf"]
    λpt = res["dual"]["ohm_pt"]
    λqt = res["dual"]["ohm_qt"]

    ωf = res["dual"]["jabr"][:, 1]
    ωt = res["dual"]["jabr"][:, 2]
    ωr = res["dual"]["jabr"][:, 3]
    ωi = res["dual"]["jabr"][:, 4]

    μθ_lb = max.(0, res["dual"]["va_diff"])
    μθ_ub = min.(0, res["dual"]["va_diff"])

    μ_w = res["dual"]["w"]
    μ_wr = res["dual"]["wr"]
    μ_wi = res["dual"]["wi"]

    # Check dual constraint corresponding to `w` variables
    δw = [
        (
            -gs[i] * λp[i]
            + bs[i] * λq[i]
            + sum(
                (gff[e] * λpf[e] - bff[e] * λqf[e] + ωf[e] / sqrt(2))
                for e in br_out[i];
                init=0
            )
            + sum(
                (gtt[e] * λpt[e] - btt[e] * λqt[e] + ωt[e] / sqrt(2))
                for e in br_in[i];
                init=0
            )
            + μ_w[i]
        )
        for i in 1:N
    ]
    @test norm(δw, Inf) <= atol

    # Check dual constraint corresponding to `wr` variables
    δwr = [
        (
            gft[e] * λpf[e] 
            + gtf[e] * λpt[e]
            - bft[e] * λqf[e] 
            - btf[e] * λqt[e]
            - tan(δθmin[e]) * μθ_lb[e]
            + tan(δθmax[e]) * μθ_ub[e]
            + ωr[e]
            + μ_wr[e]
        )
        for e in 1:E
    ]
    @test norm(δwr, Inf) <= atol

    # Check dual constraint corresponding to `wi` variables
    δwi = [
        (
            bft[e] * λpf[e] 
            - btf[e] * λpt[e]
            + gft[e] * λqf[e] 
            - gtf[e] * λqt[e]
            + μθ_lb[e]
            - μθ_ub[e]
            + ωi[e]
            + μ_wi[e]
        )
        for e in 1:E
    ]
    @test norm(δwi, Inf) <= atol
    return nothing
end

function _test_socwr_DualSolFormat()
    data = make_basic_network(pglib("pglib_opf_case118_ieee"))
    N = length(data["bus"])
    E = length(data["branch"])

    solver = CLRBL_SOLVER
    opf = PGLearn.build_opf(PGLearn.SOCOPF, data, solver)
    set_silent(opf.model)
    PGLearn.solve!(opf)

    # Check shape of dual solution
    res = PGLearn.extract_result(opf)

    @test Set(collect(keys(res))) == Set(["meta", "primal", "dual"])
    @test size(res["dual"]["jabr"]) == (E, 4)
    @test size(res["dual"]["sm_fr"]) == (E, 3)
    @test size(res["dual"]["sm_to"]) == (E, 3)
    return nothing
end

function _test_socwr128(data::Dict)
    opf = PGLearn.build_opf(PGLearn.SOCOPF, data, CLRBL128_SOLVER; T=Float128)

    PGLearn.solve!(opf)

    res = PGLearn.extract_result(opf)

    return nothing
end
