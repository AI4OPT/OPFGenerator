struct DCOPF <: AbstractFormulation end

"""
    build_opf(DCOPF, data, optimizer)

Build a DCOPF model.
"""
function build_opf(::Type{DCOPF}, data::OPFData, optimizer;
    T=Float64,    
)
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
    @variable(model, pg[1:G])

    # Branch flows
    @variable(model, pf[1:E])

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
        kcl_p[i in 1:N],
        sum(gen_status[g] * pg[g] for g in bus_gens[i])
        - sum(branch_status[a] * pf[a] for a in bus_arcs_fr[i])
        - sum(branch_status[a] * pt[a] for a in bus_arcs_to[i])  # pt == -pf
        == 
        sum(pd[l] for l in data.bus_loads[i]) + gs[i]
    )

    @constraint(model, 
        ohm_pf[e in 1:E],
        branch_status[e] * (
            -b[e] * (va[bus_fr[e]] - va[bus_to[e]])
        ) - pf[e] == 0
    )
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

    return OPFModel{DCOPF}(data, model)
end

function extract_primal(opf::OPFModel{DCOPF})
    model = opf.model
    T = JuMP.value_type(typeof(model))

    data = opf.data

    N, E, G = data.N, data.E, data.G

    primal_solution = Dict{String,Any}(
        # bus
        "va" => zeros(T, N),
        # generator
        "pg" => zeros(T, G),
        # branch
        "pf" => zeros(T, E),
    )
    if has_values(model)
        # bus
        primal_solution["va"] = value.(model[:va])
        # generator
        primal_solution["pg"] = value.(model[:pg])
        # branch
        primal_solution["pf"] = value.(model[:pf])
    end

    return primal_solution
end

function extract_dual(opf::OPFModel{DCOPF})
    model = opf.model
    T = JuMP.value_type(typeof(model))

    data = opf.data

    N, E, G = data.N, data.E, data.G

    dual_solution = Dict{String,Any}(
        # global
        "slack_bus" => zero(T),
        # bus
        "kcl_p"     => zeros(T, N),
        # generator
        # N/A
        # branch
        "ohm_pf"    => zeros(T, E),
        "va_diff"   => zeros(T, E),
        # Variable lower/upper bound
        # bus
        # N/A
        # generator
        "pg" => zeros(T, G),
        # branch
        "pf" => zeros(T, E),
    )
    if has_duals(model)
        dual_solution["slack_bus"] = dual(model[:slack_bus])

        # Bus-level constraints
        dual_solution["kcl_p"] = dual.(model[:kcl_p])

        # Generator-level constraints
        # N/A

        # Branch-level constraints
        dual_solution["ohm_pf"] = dual.(model[:ohm_pf])
        dual_solution["va_diff"] = dual.(model[:va_diff])

        # Duals of variable lower/upper bounds
        # We store λ = λₗ + λᵤ, where λₗ, λᵤ are the dual variables associated to
        #   lower and upper bounds, respectively.
        # Recall that, in JuMP's convention, we have λₗ ≥ 0, λᵤ ≤ 0, hence
        #   λₗ = max(λ, 0) and λᵤ = min(λ, 0).

        # (no bounded bus-level variables)
        # generator
        dual_solution["pg"] .+= dual.(LowerBoundRef.(model[:pg]))
        dual_solution["pg"] .+= dual.(UpperBoundRef.(model[:pg]))
        # branch
        dual_solution["pf"] .+= dual.(LowerBoundRef.(model[:pf]))
        dual_solution["pf"] .+= dual.(UpperBoundRef.(model[:pf]))
    end

    return dual_solution
end
