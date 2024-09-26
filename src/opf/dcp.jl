struct DCOPF <: AbstractFormulation end

"""
    build_opf(DCOPF, data, optimizer)

Build a DCOPF model.
"""
function build_opf(::Type{DCOPF}, network::Dict{String,Any}, optimizer;
    T=Float64,    
)
    # TODO: remove when all formulations are done
    data = OPFData(network)

    # Grab some data
    N, E, G = data.N, data.E, data.G
    gs = data.gs
    pd = data.pd
    i0 = data.ref_bus
    bus_arcs_fr, bus_arcs_to = data.bus_arcs_fr, data.bus_arcs_to
    bus_gens = data.bus_gens
    pgmin, pgmax = data.pgmin, data.pgmax
    c0, c1, c2 = data.c0, data.c1, data.c2
    gen_status = data.gen_status
    bus_fr, bus_to, b = data.bus_fr, data.bus_to, data.b
    dvamin, dvamax, smax = data.dvamin, data.dvamax, data.smax
    branch_status = data.branch_status

    model = JuMP.GenericModel{T}(optimizer)
    model.ext[:opf_model] = DCOPF

    #
    #   I. Variables
    #
    
    # nodal voltage
    @variable(model, va[1:N])

    # Active and reactive dispatch
    @variable(model, pg[g in 1:G])

    # Branch flows
    @variable(model, pf[e in 1:E])

    # 
    #   II. Constraints
    #

    # Generation bounds (both zero if generator is off)
    set_lower_bound.(pg, gen_status .* pgmin)
    set_upper_bound.(pg, gen_status .* pgmax)

    # Flow bounds (both zero if branch is off)
    set_lower_bound.(pf, branch_status .* -smax)
    set_upper_bound.(pf, branch_status .* smax)

    # Slack bus
    @constraint(model, slack_bus, va[i0] == 0.0)

    # Nodal power balance
    @expression(model, pt[e in 1:E], -pf[e])
    @constraint(model,
        kirchhoff[i in 1:N],
        sum(pg[g] for g in bus_gens[i] if gen_status[g])
        - sum(pf[a] for a in bus_arcs_fr[i] if branch_status[a])
        - sum(pt[a] for a in bus_arcs_to[i] if branch_status[a])
        == 
        pd[i] + gs[i]
    )

    model[:ohm] = Vector{ConstraintRef}(undef, E)
    model[:va_diff] = Vector{ConstraintRef}(undef, E)

    for e in 1:E if branch_status[e]
            # Branch power flow physics and limit constraints
            model[:ohm][e] = @constraint(model, -b[e] * (va[bus_fr[e]] - va[bus_to[e]]) - pf[e] == 0)

            # Voltage angle difference limit
            model[:va_diff][e] = @constraint(model, dvamin[e] ≤ va[bus_fr[e]] - va[bus_to[e]] ≤ dvamax[e])
        end
    end

    #
    #   III. Objective
    #
    l, u = extrema(c2[g] for g in 1:G if gen_status[g])
    (l == u == 0.0) || @warn "Data $(data.case) has quadratic cost terms; those terms are being ignored"

    @objective(model,
        Min,
        sum(c1[g] * pg[g] + c0[g] for g in 1:G if gen_status[g])
    )

    # TODO: update to store OPFData when all formulations are done
    return OPFModel{DCOPF}(network, model)
end

function extract_primal(opf::OPFModel{DCOPF})
    model = opf.model

    # TODO: remove when all formulations are done
    network = opf.data
    data = OPFData(network)

    N, E, G = data.N, data.E, data.G

    primal_solution = Dict{String,Any}(
        "va" => zeros(Float64, N),
        "pg" => zeros(Float64, G),
        "pf" => zeros(Float64, E),
    )
    if has_values(model)
        for i in 1:N
            primal_solution["va"][i] = value(model[:va][i])
        end

        for g in 1:G if data.gen_status[g]
                primal_solution["pg"][g] = value(model[:pg][g])
            end
        end

        for e in 1:E if data.branch_status[e]
                primal_solution["pf"][e] = value(model[:pf][e])
            end
        end
    end

    return primal_solution
end

function extract_dual(opf::OPFModel{DCOPF})
    model = opf.model

    # TODO: remove when all formulations are done
    network = opf.data
    data = OPFData(network)

    N, E, G = data.N, data.E, data.G

    dual_solution = Dict{String,Any}(
        "kirchhoff"          => zeros(Float64, N),
        "pg_lb"              => zeros(Float64, G),
        "pg_ub"              => zeros(Float64, G),
        "pf_lb"              => zeros(Float64, E),
        "pf_ub"              => zeros(Float64, E),
        "ohm"                => zeros(Float64, E),
        "va_diff"            => zeros(Float64, E),
        "slack_bus"          => 0.0,
    )
    if has_duals(model)
        for i in 1:N
            dual_solution["kirchhoff"][i] = dual(model[:kirchhoff][i])
        end

        for g in 1:G if data.gen_status[g]
                dual_solution["pg_lb"][g] = dual(LowerBoundRef(model[:pg][g]))
                dual_solution["pg_ub"][g] = dual(UpperBoundRef(model[:pg][g]))
            end
        end

        for e in 1:E if data.branch_status[e]
                dual_solution["pf_lb"][e] = dual(LowerBoundRef(model[:pf][e]))
                dual_solution["pf_ub"][e] = dual(UpperBoundRef(model[:pf][e]))
                dual_solution["ohm"][e] = dual(model[:ohm][e])
                dual_solution["va_diff"][e] = dual(model[:va_diff][e])
            end
        end

        dual_solution["slack_bus"] = dual(model[:slack_bus])
    end

    return dual_solution
end