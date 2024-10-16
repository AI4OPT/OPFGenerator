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
    @variable(model, vm[i in 1:N], start=1.0)
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
        kirchhoff_active[i in 1:N],
        sum(gen_status[g] * pg[g] for g in bus_gens[i])
        - sum(branch_status[e] * pf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * pt[e] for e in bus_arcs_to[i])
        - gs[i] * vm[i]^2
        == 
        pd[i]
    )
    @constraint(model,
        kirchhoff_reactive[i in 1:N],
        sum(gen_status[g] * qg[g] for g in bus_gens[i])
        - sum(branch_status[e] * qf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * qt[e] for e in bus_arcs_to[i])
        + bs[i] * vm[i]^2
        ==
        qd[i]
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

    @objective(model,
        Min,
        sum(c1[g] * pg[g] + c0[g] for g in 1:G if gen_status[g])
    )

    # TODO: update to store OPFData when all formulations are done
    return OPFModel{ACOPF}(network, model)
end

function extract_primal(opf::OPFModel{ACOPF})
    model = opf.model

    # TODO: remove when all formulations are done
    network = opf.data
    data = OPFData(network)

    N, E, G = data.N, data.E, data.G

    primal_solution = Dict{String,Any}(
        "vm" => zeros(Float64, N),
        "va" => zeros(Float64, N),
        "pg" => zeros(Float64, G),
        "qg" => zeros(Float64, G),
        "pf" => zeros(Float64, E),
        "qf" => zeros(Float64, E),
        "pt" => zeros(Float64, E),
        "qt" => zeros(Float64, E),
    )
    if has_values(model)
        for i in 1:N
            primal_solution["vm"][i] = value(model[:vm][i])
            primal_solution["va"][i] = value(model[:va][i])
        end

        for g in 1:G if data.gen_status[g]
                primal_solution["pg"][g] = value(model[:pg][g])
                primal_solution["qg"][g] = value(model[:qg][g])
            end
        end

        for e in 1:E if data.branch_status[e]
                primal_solution["pf"][e] = value(model[:pf][e])
                primal_solution["qf"][e] = value(model[:qf][e])
                primal_solution["pt"][e] = value(model[:pt][e])
                primal_solution["qt"][e] = value(model[:qt][e])
            end
        end
    end

    return primal_solution
end

function extract_dual(opf::OPFModel{ACOPF})
    model = opf.model

    # TODO: remove when all formulations are done
    network = opf.data
    data = OPFData(network)

    N, E, G = data.N, data.E, data.G

    dual_solution = Dict{String,Any}(
        "vm_lb"              => zeros(Float64, N),
        "vm_ub"              => zeros(Float64, N),
        "kirchhoff_active"   => zeros(Float64, N),
        "kirchhoff_reactive" => zeros(Float64, N),
        "pg_lb"              => zeros(Float64, G),
        "pg_ub"              => zeros(Float64, G),
        "qg_lb"              => zeros(Float64, G),
        "qg_ub"              => zeros(Float64, G),
        "sm_fr"              => zeros(Float64, E),
        "sm_to"              => zeros(Float64, E),
        "ohm_pf"             => zeros(Float64, E),
        "ohm_pt"             => zeros(Float64, E),
        "ohm_qf"             => zeros(Float64, E),
        "ohm_qt"             => zeros(Float64, E),
        "va_diff"            => zeros(Float64, E),
        "slack_bus"          => 0.0,
    )
    if has_duals(model)
        for i in 1:N
            dual_solution["vm_lb"][i] = dual(LowerBoundRef(model[:vm][i]))
            dual_solution["vm_ub"][i] = dual(UpperBoundRef(model[:vm][i]))
            dual_solution["kirchhoff_active"][i] = dual(model[:kirchhoff_active][i])
            dual_solution["kirchhoff_reactive"][i] = dual(model[:kirchhoff_reactive][i])
        end

        for g in 1:G if data.gen_status[g]
                dual_solution["pg_lb"][g] = dual(LowerBoundRef(model[:pg][g]))
                dual_solution["pg_ub"][g] = dual(UpperBoundRef(model[:pg][g]))
                dual_solution["qg_lb"][g] = dual(LowerBoundRef(model[:qg][g]))
                dual_solution["qg_ub"][g] = dual(UpperBoundRef(model[:qg][g]))
            end
        end

        for e in 1:E if data.branch_status[e]
                dual_solution["sm_fr"][e] = dual(model[:sm_fr][e])
                dual_solution["sm_to"][e] = dual(model[:sm_to][e])
                dual_solution["ohm_pf"][e] = dual(model[:ohm_pf][e])
                dual_solution["ohm_pt"][e] = dual(model[:ohm_pt][e])
                dual_solution["ohm_qf"][e] = dual(model[:ohm_qf][e])
                dual_solution["ohm_qt"][e] = dual(model[:ohm_qt][e])
                dual_solution["va_diff"][e] = dual(model[:va_diff][e])
            end
        end

        dual_solution["slack_bus"] = dual(model[:slack_bus])
    end

    return dual_solution
end