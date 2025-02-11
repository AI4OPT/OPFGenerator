struct SDPOPF <: AbstractFormulation end

"""
    build_opf(SDPOPF, data, optimizer)

Build an SDPOPF model.
"""
function build_opf(::Type{SDPOPF}, data::OPFData, optimizer;
    T=Float64,
)
    # Grab some data
    N, E, G = data.N, data.E, data.G
    vmin, vmax = data.vmin, data.vmax
    i0 = data.ref_bus
    gs, bs = data.gs, data.bs
    pd, qd = data.pd, data.qd
    bus_arcs_fr, bus_arcs_to = data.bus_arcs_fr, data.bus_arcs_to
    bus_gens = data.bus_gens
    pgmin, pgmax = data.pgmin, data.pgmax
    qgmin, qgmax = data.qgmin, data.qgmax
    c0, c1, c2 = data.c0, data.c1, data.c2
    gen_status = data.gen_status
    bus_fr, bus_to = data.bus_fr, data.bus_to
    gff, gft, gtf, gtt = data.gff, data.gft, data.gtf, data.gtt
    bff, bft, btf, btt = data.bff, data.bft, data.btf, data.btt
    dvamin, dvamax, smax = data.dvamin, data.dvamax, data.smax
    branch_status = data.branch_status

    wr_min, wr_max, wi_min, wi_max = compute_voltage_phasor_bounds(data)

    model = JuMP.GenericModel{T}(optimizer)
    model.ext[:opf_model] = SDPOPF

    #
    #   I. Variables
    #
    #=
        Some generators and branches may be inactive, as indicated by `branch_status` and `gen_status`.
        Primal variables associated to inactive components are still declared, and handled as follows:
        * lower/upper bounds are set to zero
        * constraint coefficients are set to zero
    =#

    # voltage magnitude and product
    WR_start = zeros(N, N) + I
    @variable(model, WR[i=1:N, j=1:N], Symmetric, start=WR_start[i,j])
    @variable(model, WI[1:N, 1:N] in SkewSymmetricMatrixSpace(), start=0.0)

    # Active and reactive dispatch
    @variable(model, pg[g in 1:G])
    @variable(model, qg[g in 1:G])
    
    # Directional branch flows
    @variable(model, pf[e in 1:E])
    @variable(model, qf[e in 1:E])
    @variable(model, pt[e in 1:E])
    @variable(model, qt[e in 1:E])

    # 
    #   II. Constraints
    #

    # Voltage magnitude bounds
    set_lower_bound.([WR[i, i] for i in 1:N], vmin.^2)
    set_upper_bound.([WR[i, i] for i in 1:N], vmax.^2)

    # Voltage product bounds
    set_lower_bound.([WR[bus_fr[e], bus_to[e]] for e in 1:E], wr_min)
    set_upper_bound.([WR[bus_fr[e], bus_to[e]] for e in 1:E], wr_max)
    set_lower_bound.([WI[sort([bus_fr[e], bus_to[e]])...] for e in 1:E], wi_min)
    set_upper_bound.([WI[sort([bus_fr[e], bus_to[e]])...] for e in 1:E], wi_max)

    # Active generation bounds (both zero if generator is off)
    set_lower_bound.(pg, gen_status .* pgmin)
    set_upper_bound.(pg, gen_status .* pgmax)

    # Reactive generation bounds (both zero if generator is off)
    set_lower_bound.(qg, gen_status .* qgmin)
    set_upper_bound.(qg, gen_status .* qgmax)

    # Active flow bounds (both zero if branch is off)
    set_lower_bound.(pf, branch_status .* -smax)
    set_upper_bound.(pf, branch_status .* smax)
    set_lower_bound.(pt, branch_status .* -smax)
    set_upper_bound.(pt, branch_status .* smax)

    # Reactive flow bounds (both zero if branch is off)
    set_lower_bound.(qf, branch_status .* -smax)
    set_upper_bound.(qf, branch_status .* smax)
    set_lower_bound.(qt, branch_status .* -smax)
    set_upper_bound.(qt, branch_status .* smax)

    # Nodal power balance
    @constraint(model, 
        kcl_p[i in 1:N],
        sum(gen_status[g] * pg[g] for g in bus_gens[i])
        - sum(branch_status[e] * pf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * pt[e] for e in bus_arcs_to[i])
        - gs[i] * WR[i, i]
        ==
        sum(pd[l] for l in data.bus_loads[i])
    )
    @constraint(model,
        kcl_q[i in 1:N],
        sum(gen_status[g] * qg[g] for g in bus_gens[i]) 
        - sum(branch_status[e] * qf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * qt[e] for e in bus_arcs_to[i])
        + bs[i] * WR[i, i]
        ==
        sum(qd[l] for l in data.bus_loads[i])
    )

    # Branch power flow physics and limit constraints
    # Note that the same variable may appear more than once if there are parallel branches
    @expression(model, wf[e in 1:E], WR[bus_fr[e], bus_fr[e]])
    @expression(model, wt[e in 1:E], WR[bus_to[e], bus_to[e]])
    @expression(model, wr[e in 1:E], WR[bus_fr[e], bus_to[e]])
    @expression(model, wi[e in 1:E], WI[bus_fr[e], bus_to[e]])

    # Ohm's law
    @constraint(model, ohm_pf[e in 1:E],
        branch_status[e] * ( gff[e] * wf[e] + gft[e] * wr[e] + bft[e] * wi[e]) - pf[e] == 0
    )
    @constraint(model, ohm_qf[e in 1:E],
        branch_status[e] * (-bff[e] * wf[e] - bft[e] * wr[e] + gft[e] * wi[e]) - qf[e] == 0
    )
    @constraint(model, ohm_pt[e in 1:E],
        branch_status[e] * ( gtt[e] * wt[e] + gtf[e] * wr[e] - btf[e] * wi[e]) - pt[e] == 0
    )
    @constraint(model, ohm_qt[e in 1:E],
        branch_status[e] * (-btt[e] * wt[e] - btf[e] * wr[e] - gtf[e] * wi[e]) - qt[e] == 0
    )

    # Thermal limit
    @constraint(model, sm_fr[e in 1:E], [smax[e], pf[e], qf[e]] in SecondOrderCone())
    @constraint(model, sm_to[e in 1:E], [smax[e], pt[e], qt[e]] in SecondOrderCone())

    # Voltage angle difference limit
    @constraint(model, va_diff_lb[e in 1:E], branch_status[e] * wi[e] - branch_status[e] * tan(dvamin[e]) * wr[e] >= 0)
    @constraint(model, va_diff_ub[e in 1:E], branch_status[e] * wi[e] - branch_status[e] * tan(dvamax[e]) * wr[e] <= 0)
    
    # PSD constraint
    @constraint(model, S, [WR WI; -WI WR] in PSDCone())
    
    #
    #   III. Objective
    #
    l, u = extrema(c2)
    (l == u == 0.0) || @warn "Data $(data.case) has quadratic cost terms; those terms are being ignored"
    @objective(model,
        Min,
        sum(c1[g] * pg[g] + c0[g] for g in 1:G if gen_status[g])
    )

    return OPFModel{SDPOPF}(data, model)
end

function extract_primal(opf::OPFModel{SDPOPF})
    model = opf.model
    T = JuMP.value_type(typeof(model))

    data = opf.data

    N, E, G = data.N, data.E, data.G

    primal_solution = Dict{String,Any}(
        # bus
        "w" => zeros(T, N),
        # generator
        "pg" => zeros(T, G),
        "qg" => zeros(T, G),
        # branch
        "wr" => zeros(T, E),
        "wi" => zeros(T, E),
        "pf" => zeros(T, E),
        "qf" => zeros(T, E),
        "pt" => zeros(T, E),
        "qt" => zeros(T, E),
    )
    if has_values(model)
        # bus
        primal_solution["w"] = value.([model[:WR][i, i] for i in 1:N])

        # generator
        primal_solution["pg"] = value.(model[:pg])
        primal_solution["qg"] = value.(model[:qg])

        # branch
        # By construction of the PSD matrix, WR and WI are defined for buspairs instead of branches.
        # For convenience, they are stored for branches.
        # The same values will be repeated if there are multiple branches between a bus pair.
        # Only extract WR and WI entries that correspond to branches (necessary for power flow computation)
        # Other entries are not extracted as they are dense
        primal_solution["wr"] = value.(model[:wr])
        primal_solution["wi"] = value.(model[:wi])
        primal_solution["pf"] = value.(model[:pf])
        primal_solution["qf"] = value.(model[:qf])
        primal_solution["pt"] = value.(model[:pt])
        primal_solution["qt"] = value.(model[:qt])
    end

    return primal_solution
end

function extract_dual(opf::OPFModel{SDPOPF})
    model = opf.model
    T = JuMP.value_type(typeof(model))

    data = opf.data

    N, E, G = data.N, data.E, data.G
    bus_fr, bus_to = data.bus_fr, data.bus_to

    dual_solution = Dict{String,Any}(
        # bus
        "kcl_p"      => zeros(T, N),
        "kcl_q"      => zeros(T, N),
        "s"          => zeros(T, N),
        # generator
        # N/A
        # branch
        "ohm_pf"     => zeros(T, E),
        "ohm_pt"     => zeros(T, E),
        "ohm_qf"     => zeros(T, E),
        "ohm_qt"     => zeros(T, E),
        "va_diff"    => zeros(T, E),
        "sm_fr"      => zeros(T, E, 3),
        "sm_to"      => zeros(T, E, 3),
        "sr"         => zeros(T, E),
        "si"         => zeros(T, E),
        # variables lower/upper bounds
        # bus
        "w"          => zeros(T, N),
        # generator
        "pg"         => zeros(T, G),
        "qg"         => zeros(T, G),
        # branch
        "wr"         => zeros(T, E),
        "wi"         => zeros(T, E),
        "pf"         => zeros(T, E),
        "qf"         => zeros(T, E),
        "pt"         => zeros(T, E),
        "qt"         => zeros(T, E),
    )

    if has_duals(model)
        S = dual.(model[:S])

        # Bus-level constraints
        dual_solution["kcl_p"] = dual.(model[:kcl_p])
        dual_solution["kcl_q"] = dual.(model[:kcl_q])
        dual_solution["s"] = [S[i, i] for i in 1:N] # upper-left diagonal of S

        # Generator-level constraints
        # N/A

        # Branch-level constraints
        dual_solution["ohm_pf"] = dual.(model[:ohm_pf])
        dual_solution["ohm_pt"] = dual.(model[:ohm_pt])
        dual_solution["ohm_qf"] = dual.(model[:ohm_qf])
        dual_solution["ohm_qt"] = dual.(model[:ohm_qt])
        dual_solution["va_diff"] = dual.(model[:va_diff_lb]) + dual.(model[:va_diff_ub])  # same as bound constraints
        dual_solution["sm_fr"] = dual.(model[:sm_fr])
        dual_solution["sm_to"] = dual.(model[:sm_to])
        # By construction of the PSD matrix, sr and si are defined for buspairs instead of branches.
        # For convenience, they are stored for branches.
        # The same values will be repeated if there are multiple branches between a bus pair.
        dual_solution["sr"] = [S[bus_fr[e], bus_to[e]] for e in 1:E] # upper left block
        dual_solution["si"] = [S[bus_fr[e], bus_to[e] + N] for e in 1:E] # upper right block
        
        # For conic constraints, JuMP will return Vector{Vector{T}}
        # reshape duals of conic constraints into matrix shape
        dual_solution["sm_fr"] = mapreduce(permutedims, vcat, dual_solution["sm_fr"])
        dual_solution["sm_to"] = mapreduce(permutedims, vcat, dual_solution["sm_to"])

        # Duals of variable lower/upper bounds
        # We store λ = λₗ + λᵤ, where λₗ, λᵤ are the dual variables associated to
        #   lower and upper bounds, respectively.
        # Recall that, in JuMP's convention, we have λₗ ≥ 0, λᵤ ≤ 0, hence
        #   λₗ = max(λ, 0) and λᵤ = min(λ, 0).

        # bus
        dual_solution["w"]  = dual.(LowerBoundRef.([model[:WR][i, i] for i in 1:N])) + dual.(UpperBoundRef.([model[:WR][i, i] for i in 1:N]))
        # generator
        dual_solution["pg"] = dual.(LowerBoundRef.(model[:pg])) + dual.(UpperBoundRef.(model[:pg]))
        dual_solution["qg"] = dual.(LowerBoundRef.(model[:qg])) + dual.(UpperBoundRef.(model[:qg]))
        # branch
        # Note that wr and wi are by construction defined on buspairs, not branches.
        # The same values will be repeated if there are multiple branches between a bus pair.
        dual_solution["wr"] = dual.(LowerBoundRef.([model[:WR][bus_fr[e], bus_to[e]] for e in 1:E])) + dual.(UpperBoundRef.([model[:WR][bus_fr[e], bus_to[e]] for e in 1:E]))
        # In the current implementation of SkewSymmetricMatrixSpace, entries of WI are AffExpr, so need to extract the corresponding variable
        dual_solution["wi"] = dual.(LowerBoundRef.([first(keys(model[:WI][bus_fr[e], bus_to[e]].terms)) for e in 1:E])) + dual.(UpperBoundRef.([first(keys(model[:WI][bus_fr[e], bus_to[e]].terms)) for e in 1:E]))
        dual_solution["pf"] = dual.(LowerBoundRef.(model[:pf])) + dual.(UpperBoundRef.(model[:pf]))
        dual_solution["qf"] = dual.(LowerBoundRef.(model[:qf])) + dual.(UpperBoundRef.(model[:qf]))
        dual_solution["pt"] = dual.(LowerBoundRef.(model[:pt])) + dual.(UpperBoundRef.(model[:pt]))
        dual_solution["qt"] = dual.(LowerBoundRef.(model[:qt])) + dual.(UpperBoundRef.(model[:qt]))
    end

    return dual_solution
end