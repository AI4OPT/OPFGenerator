struct ACOPF <: AbstractFormulation end

"""
    build_opf(ACOPF, data, optimizer)

Build an ACOPF model.
"""
function build_opf(::Type{ACOPF}, network::Dict{String,Any}, optimizer;
    T=Float64,    
)
    # TODO: remove when all formulations are done
    data = OPFData(network)
    
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

    model = JuMP.GenericModel{T}(optimizer)
    model.ext[:opf_model] = ACOPF

    #
    #   I. Variables
    #
    
    # nodal voltage
    @variable(model, vm[1:N], start=1.0)
    @variable(model, va[1:N])

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
    set_lower_bound.(vm, vmin)
    set_upper_bound.(vm, vmax)

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

    # Slack bus
    @constraint(model, slack_bus, va[i0] == 0.0)

    # Nodal power balance
    @constraint(model,
        kcl_p[i in 1:N],
        sum(gen_status[g] * pg[g] for g in bus_gens[i])
        - sum(branch_status[e] * pf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * pt[e] for e in bus_arcs_to[i])
        - gs[i] * vm[i]^2
        == 
        sum(pd[l] for l in data.bus_loads[i])
    )
    @constraint(model,
        kcl_q[i in 1:N],
        sum(gen_status[g] * qg[g] for g in bus_gens[i])
        - sum(branch_status[e] * qf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * qt[e] for e in bus_arcs_to[i])
        + bs[i] * vm[i]^2
        ==
        sum(qd[l] for l in data.bus_loads[i])
    )


    # Ohm's law
    # Some useful expressions first
    @expression(model, wf[e in 1:E], vm[bus_fr[e]]^2)
    @expression(model, wt[e in 1:E], vm[bus_to[e]]^2)
    @expression(model, wr[e in 1:E], vm[bus_fr[e]] * vm[bus_to[e]] * cos(va[bus_fr[e]] - va[bus_to[e]]))
    @expression(model, wi[e in 1:E], vm[bus_fr[e]] * vm[bus_to[e]] * sin(va[bus_fr[e]] - va[bus_to[e]]))
    # Actual constraints
    @constraint(model,
        ohm_pf[e in 1:E],
        branch_status[e] * ( gff[e] * wf[e] + gft[e] * wr[e] + bft[e] * wi[e]) - pf[e] == 0
    )
    @constraint(model,
        ohm_qf[e in 1:E],
        branch_status[e] * (-bff[e] * wf[e] - bft[e] * wr[e] + gft[e] * wi[e]) - qf[e] == 0
    )
    @constraint(model,
        ohm_pt[e in 1:E],
        branch_status[e] * ( gtt[e] * wt[e] + gtf[e] * wr[e] - btf[e] * wi[e]) - pt[e] == 0
    )
    @constraint(model,
        ohm_qt[e in 1:E],
        branch_status[e] * (-btt[e] * wt[e] - btf[e] * wr[e] - gtf[e] * wi[e]) - qt[e] == 0
    )
    
    # Thermal limit
    @constraint(model, sm_fr[e in 1:E], pf[e]^2 + qf[e]^2 ≤ smax[e]^2)
    @constraint(model, sm_to[e in 1:E], pt[e]^2 + qt[e]^2 ≤ smax[e]^2)

    # Voltage angle difference limit
    @constraint(model,
        va_diff[e in 1:E],
        dvamin[e] ≤ branch_status[e] * (va[bus_fr[e]] - va[bus_to[e]]) ≤ dvamax[e]
    )

    #
    #   III. Objective
    #
    l, u = extrema(c2)
    (l == u == 0.0) || @warn "Data $(data.case) has quadratic cost terms; those terms are being ignored"
    @objective(model,
        Min,
        sum(c1[g] * pg[g] + c0[g] for g in 1:G if gen_status[g])
    )

    # TODO: update to store OPFData when all formulations are done
    return OPFModel{ACOPF}(network, model)
end

function extract_primal(opf::OPFModel{ACOPF})
    model = opf.model
    T = JuMP.value_type(typeof(model))

    # TODO: remove when all formulations are done
    network = opf.data
    data = OPFData(network)

    N, E, G = data.N, data.E, data.G

    primal_solution = Dict{String,Any}(
        # bus
        "vm" => zeros(T, N),
        "va" => zeros(T, N),
        # generator
        "pg" => zeros(T, G),
        "qg" => zeros(T, G),
        # branch
        "pf" => zeros(T, E),
        "qf" => zeros(T, E),
        "pt" => zeros(T, E),
        "qt" => zeros(T, E),
    )
    if has_values(model)
        # bus
        primal_solution["vm"] = value.(model[:vm])
        primal_solution["va"] = value.(model[:va])

        # generator
        primal_solution["pg"] = value.(model[:pg])
        primal_solution["qg"] = value.(model[:qg])

        # branch
        primal_solution["pf"] = value.(model[:pf])
        primal_solution["qf"] = value.(model[:qf])
        primal_solution["pt"] = value.(model[:pt])
        primal_solution["qt"] = value.(model[:qt])
    end

    return primal_solution
end

function extract_dual(opf::OPFModel{ACOPF})
    model = opf.model
    T = JuMP.value_type(typeof(model))

    # TODO: remove when all formulations are done
    network = opf.data
    data = OPFData(network)

    N, E, G = data.N, data.E, data.G

    dual_solution = Dict{String,Any}(
        # global
        "slack_bus" => zero(T),
        # bus
        "kcl_p"     => zeros(T, N),
        "kcl_q"     => zeros(T, N),
        # generator
        # N/A
        # branch
        "ohm_pf"    => zeros(T, E),
        "ohm_pt"    => zeros(T, E),
        "ohm_qf"    => zeros(T, E),
        "ohm_qt"    => zeros(T, E),
        "sm_fr"     => zeros(T, E),
        "sm_to"     => zeros(T, E),
        "va_diff"   => zeros(T, E),
        # variables lower/upper bounds
        # bus
        "vm_lb"     => zeros(T, N),
        "vm_ub"     => zeros(T, N),
        # generator
        "pg_lb"     => zeros(T, G),
        "pg_ub"     => zeros(T, G),
        "qg_lb"     => zeros(T, G),
        "qg_ub"     => zeros(T, G),
        # branch
        "pf_lb"      => zeros(T, E),
        "pf_ub"      => zeros(T, E),
        "qf_lb"      => zeros(T, E),
        "qf_ub"      => zeros(T, E),
        "pt_lb"      => zeros(T, E),
        "pt_ub"      => zeros(T, E),
        "qt_lb"      => zeros(T, E),
        "qt_ub"      => zeros(T, E),
    )

    if has_duals(model)
        # global
        dual_solution["slack_bus"] = dual(model[:slack_bus])
        # bus
        dual_solution["kcl_p"] = dual.(model[:kcl_p])
        dual_solution["kcl_q"] = dual.(model[:kcl_q])

        # generator
        # N/A

        # branch
        dual_solution["ohm_pf"] = dual.(model[:ohm_pf])
        dual_solution["ohm_pt"] = dual.(model[:ohm_pt])
        dual_solution["ohm_qf"] = dual.(model[:ohm_qf])
        dual_solution["ohm_qt"] = dual.(model[:ohm_qt])
        dual_solution["sm_fr"] = dual.(model[:sm_fr])
        dual_solution["sm_to"] = dual.(model[:sm_to])
        dual_solution["va_diff"] = dual.(model[:va_diff])

        # Variable lower/upper bounds
        # bus
        dual_solution["vm_lb"] = dual.(LowerBoundRef.(model[:vm]))
        dual_solution["vm_ub"] = dual.(UpperBoundRef.(model[:vm]))
        #    (nodal voltage angles have no lower/upper bounds)
        # generator
        dual_solution["pg_lb"] = dual.(LowerBoundRef.(model[:pg]))
        dual_solution["pg_ub"] = dual.(UpperBoundRef.(model[:pg]))
        dual_solution["qg_lb"] = dual.(LowerBoundRef.(model[:qg]))
        dual_solution["qg_ub"] = dual.(UpperBoundRef.(model[:qg]))
        # branch
        dual_solution["pf_lb"] = dual.(LowerBoundRef.(model[:pf]))
        dual_solution["pf_ub"] = dual.(UpperBoundRef.(model[:pf]))
        dual_solution["qf_lb"] = dual.(LowerBoundRef.(model[:qf]))
        dual_solution["qf_ub"] = dual.(UpperBoundRef.(model[:qf]))
        dual_solution["pt_lb"] = dual.(LowerBoundRef.(model[:pt]))
        dual_solution["pt_ub"] = dual.(UpperBoundRef.(model[:pt]))
        dual_solution["qt_lb"] = dual.(LowerBoundRef.(model[:qt]))
        dual_solution["qt_ub"] = dual.(UpperBoundRef.(model[:qt]))
    end

    return dual_solution
end
