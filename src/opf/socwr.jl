struct SOCOPFQuad <: AbstractFormulation end
struct SOCOPF <: AbstractFormulation end

"""
    build_opf(SOCOPF, data, optimizer)

Build an SOCOPF model.
"""
function build_opf(::Type{OPF}, network::Dict{String,Any}, optimizer;
    T=Float64,    
) where {OPF <: Union{SOCOPFQuad,SOCOPF}}
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

    wr_min, wr_max, wi_min, wi_max = compute_voltage_phasor_bounds(data)

    model = JuMP.GenericModel{T}(optimizer)
    model.ext[:opf_model] = OPF

    #
    #   I. Variables
    #
    
    # voltage magnitude and product
    @variable(model, w[1:N])
    @variable(model, wr[1:E])
    @variable(model, wi[1:E])

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
    set_lower_bound.(w, vmin.^2)
    set_upper_bound.(w, vmax.^2)

    # Voltage product bounds
    set_lower_bound.(wr, wr_min)
    set_upper_bound.(wr, wr_max)
    set_lower_bound.(wi, wi_min)
    set_upper_bound.(wi, wi_max)

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
        kirchhoff_active[i in 1:N],
        sum(pg[g] for g in bus_gens[i] if gen_status[g])
        - sum(pf[e] for e in bus_arcs_fr[i] if branch_status[e])
        - sum(pt[e] for e in bus_arcs_to[i] if branch_status[e])
        - gs[i] * w[i]
        == 
        pd[i]
    )
    @constraint(model,
        kirchhoff_reactive[i in 1:N],
        sum(qg[g] for g in bus_gens[i] if gen_status[g])
        - sum(qf[e] for e in bus_arcs_fr[i] if branch_status[e])
        - sum(qt[e] for e in bus_arcs_to[i] if branch_status[e])
        + bs[i] * w[i]
        ==
        qd[i]
    )

    # Branch power flow physics and limit constraints
    @expression(model, wf[e in 1:E], w[bus_fr[e]])
    @expression(model, wt[e in 1:E], w[bus_to[e]])

    model[:sm_fr] = Vector{ConstraintRef}(undef, E)
    model[:sm_to] = Vector{ConstraintRef}(undef, E)
    model[:va_diff_lb] = Vector{ConstraintRef}(undef, E)
    model[:va_diff_ub] = Vector{ConstraintRef}(undef, E)
    model[:ohm_pf] = Vector{ConstraintRef}(undef, E)
    model[:ohm_qf] = Vector{ConstraintRef}(undef, E)
    model[:ohm_pt] = Vector{ConstraintRef}(undef, E)
    model[:ohm_qt] = Vector{ConstraintRef}(undef, E)
    model[:jabr] = Vector{ConstraintRef}(undef, E)

    for e in 1:E if branch_status[e]
            # Ohm's law
            model[:ohm_pf][e] = @constraint(model,  gff[e] * wf[e] + gft[e] * wr[e] + bft[e] * wi[e] - pf[e] == 0)
            model[:ohm_qf][e] = @constraint(model, -bff[e] * wf[e] - bft[e] * wr[e] + gft[e] * wi[e] - qf[e] == 0)
            model[:ohm_pt][e] = @constraint(model,  gtt[e] * wt[e] + gtf[e] * wr[e] - btf[e] * wi[e] - pt[e] == 0)
            model[:ohm_qt][e] = @constraint(model, -btt[e] * wt[e] - btf[e] * wr[e] - gtf[e] * wi[e] - qt[e] == 0)
        
            # Thermal limit
            if OPF == SOCOPF
                model[:sm_fr][e] = @constraint(model, [smax[e], pf[e], qf[e]] in SecondOrderCone())
                model[:sm_to][e] = @constraint(model, [smax[e], pt[e], qt[e]] in SecondOrderCone())
            elseif OPF == SOCOPFQuad
                model[:sm_fr][e] = @constraint(model, pf[e]^2 + qf[e]^2 ≤ smax[e]^2)
                model[:sm_to][e] = @constraint(model, pt[e]^2 + qt[e]^2 ≤ smax[e]^2)
            end
        
            # Voltage angle difference limit
            model[:va_diff_lb][e] = @constraint(model, wi[e] >= tan(dvamin[e]) * wr[e])
            model[:va_diff_ub][e] = @constraint(model, wi[e] <= tan(dvamax[e]) * wr[e])

            # Jabr constraints
            if OPF == SOCOPF
                model[:jabr][e] = @constraint(model, [wf[e] / sqrt(2), wt[e] / sqrt(2), wr[e], wi[e]] in RotatedSecondOrderCone())
            elseif OPF == SOCOPFQuad
                model[:jabr][e] = @constraint(model, wr[e]^2 + wi[e]^2 ≤ wf[e] * wt[e])
            end
        end
    end

    #
    #   III. Objective
    #

    @objective(model,
        Min,
        sum(c2[g] * pg[g]^2 + c1[g] * pg[g] + c0[g] for g in 1:G if gen_status[g])
    )

    # TODO: update to store OPFData when all formulations are done
    return OPFModel{SOCOPF}(network, model)
end

function extract_primal(opf::OPFModel{OPF}) where {OPF <: Union{SOCOPFQuad,SOCOPF}}
    model = opf.model

    # TODO: remove when all formulations are done
    network = opf.data
    data = OPFData(network)

    N, E, G = data.N, data.E, data.G

    primal_solution = Dict{String,Any}(
        "w" => zeros(Float64, N),
        "wr" => zeros(Float64, E),
        "wi" => zeros(Float64, E),
        "pg" => zeros(Float64, G),
        "qg" => zeros(Float64, G),
        "pf" => zeros(Float64, E),
        "qf" => zeros(Float64, E),
        "pt" => zeros(Float64, E),
        "qt" => zeros(Float64, E),
    )
    if has_values(model)
        for i in 1:N
            primal_solution["w"][i] = value(model[:w][i])
        end

        for g in 1:G if data.gen_status[g]
                primal_solution["pg"][g] = value(model[:pg][g])
                primal_solution["qg"][g] = value(model[:qg][g])
            end
        end

        for e in 1:E if data.branch_status[e]
                primal_solution["wr"][e] = value(model[:wr][e])
                primal_solution["wi"][e] = value(model[:wi][e])
                primal_solution["pf"][e] = value(model[:pf][e])
                primal_solution["qf"][e] = value(model[:qf][e])
                primal_solution["pt"][e] = value(model[:pt][e])
                primal_solution["qt"][e] = value(model[:qt][e])
            end
        end
    end

    return primal_solution
end

function extract_dual(opf::OPFModel{OPF}) where {OPF <: Union{SOCOPFQuad,SOCOPF}}
    model = opf.model

    # TODO: remove when all formulations are done
    network = opf.data
    data = OPFData(network)

    N, E, G = data.N, data.E, data.G

    dual_solution = Dict{String,Any}(
        "w_lb"               => zeros(Float64, N),
        "w_ub"               => zeros(Float64, N),
        "kirchhoff_active"   => zeros(Float64, N),
        "kirchhoff_reactive" => zeros(Float64, N),
        "pg_lb"              => zeros(Float64, G),
        "pg_ub"              => zeros(Float64, G),
        "qg_lb"              => zeros(Float64, G),
        "qg_ub"              => zeros(Float64, G),
        "wr_lb"              => zeros(Float64, E),
        "wr_ub"              => zeros(Float64, E),
        "wi_lb"              => zeros(Float64, E),
        "wi_ub"              => zeros(Float64, E),
        "ohm_pf"             => zeros(Float64, E),
        "ohm_pt"             => zeros(Float64, E),
        "ohm_qf"             => zeros(Float64, E),
        "ohm_qt"             => zeros(Float64, E),
        "va_diff_lb"         => zeros(Float64, E),
        "va_diff_ub"         => zeros(Float64, E),
    )

    if OPF == SOCOPFQuad
        dual_solution["sm_fr"] = zeros(Float64, E)
        dual_solution["sm_to"] = zeros(Float64, E)
        dual_solution["jabr"] = zeros(Float64, E)
    elseif OPF == SOCOPF
        dual_solution["sm_fr"] = zeros(Float64, E, 3)
        dual_solution["sm_to"] = zeros(Float64, E, 3)
        dual_solution["jabr"] = zeros(Float64, E, 4)
    end

    if has_duals(model)
        for i in 1:N
            dual_solution["w_lb"][i] = dual(LowerBoundRef(model[:w][i]))
            dual_solution["w_ub"][i] = dual(UpperBoundRef(model[:w][i]))
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
                dual_solution["wr_lb"][e] = dual(LowerBoundRef(model[:wr][e]))
                dual_solution["wr_ub"][e] = dual(UpperBoundRef(model[:wr][e]))
                dual_solution["wi_lb"][e] = dual(LowerBoundRef(model[:wi][e]))
                dual_solution["wi_ub"][e] = dual(UpperBoundRef(model[:wi][e]))
                dual_solution["ohm_pf"][e] = dual(model[:ohm_pf][e])
                dual_solution["ohm_pt"][e] = dual(model[:ohm_pt][e])
                dual_solution["ohm_qf"][e] = dual(model[:ohm_qf][e])
                dual_solution["ohm_qt"][e] = dual(model[:ohm_qt][e])
                dual_solution["va_diff_lb"][e] = dual(model[:va_diff_lb][e])
                dual_solution["va_diff_ub"][e] = dual(model[:va_diff_ub][e])

                if OPF == SOCOPFQuad
                    dual_solution["sm_fr"][e] = dual(model[:sm_fr][e])
                    dual_solution["sm_to"][e] = dual(model[:sm_to][e])
                    dual_solution["jabr"][e] = dual(model[:jabr][e])
                elseif OPF == SOCOPF
                    dual_solution["sm_fr"][e, :] .= dual(model[:sm_fr][e])
                    dual_solution["sm_to"][e, :] .= dual(model[:sm_to][e])
                    dual_solution["jabr"][e, :] .= dual(model[:jabr][e])
                end
            end
        end
    end

    return dual_solution
end