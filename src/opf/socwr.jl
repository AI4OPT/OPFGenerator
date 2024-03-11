"""
    build_soc_opf(data, optimizer)

Build an SOC-OPF model.

This implementation is based on the SOC-OPF formulation of PMAnnex.jl
    https://github.com/lanl-ansi/PMAnnex.jl/blob/f303f3c3c61e2d1a050ee7651fa6e8abc4055b55/src/model/opf.jl
"""
function build_opf(::Type{OPF}, data::Dict{String,Any}, optimizer;
    T=Float64,
    wbounds=false,
) where {OPF <: Union{PM.SOCWRPowerModel,PM.SOCWRConicPowerModel}}
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    # Cleanup and pre-process data
    PM.standardize_cost_terms!(data, order=2)
    PM.calc_thermal_limits!(data)
    ref = PM.build_ref(data)[:it][:pm][:nw][0]
    
    # Grab some data
    N = length(ref[:bus])
    G = length(ref[:gen])
    E = length(ref[:branch])
    bus_loads = [
        [ref[:load][l] for l in ref[:bus_loads][i]]
        for i in 1:N
    ]
    bus_shunts = [
        [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        for i in 1:N
    ]
    # Branch --> bus-pair correspondence
    br2bp = [
        (ref[:branch][e]["f_bus"], ref[:branch][e]["t_bus"])
        for e in 1:E
    ]

    model = JuMP.GenericModel{T}(optimizer)
    model.ext[:opf_model] = OPF

    #
    #   I. Variables
    #
    # nodal voltage
    JuMP.@variable(model, ref[:bus][i]["vmin"]^2 <= w[i in 1:N] <= ref[:bus][i]["vmax"]^2, start=1.001)

    # wr, wi variables (per branch)
    JuMP.@variable(model, wr[e in 1:E], start=1.0)
    JuMP.@variable(model, wi[e in 1:E])
    model.ext[:wbounds] = wbounds
    wr_min, wr_max, wi_min, wi_max = PM.ref_calc_voltage_product_bounds(ref[:buspairs])
    if wbounds
        for e in 1:E
            set_lower_bound(wr[e], wr_min[br2bp[e]])
            set_upper_bound(wr[e], wr_max[br2bp[e]])
            set_lower_bound(wi[e], wi_min[br2bp[e]])
            set_upper_bound(wi[e], wi_max[br2bp[e]])
        end
    end
    # Active and reactive dispatch
    JuMP.@variable(model, ref[:gen][g]["pmin"] <= pg[g in 1:G] <= ref[:gen][g]["pmax"])
    JuMP.@variable(model, ref[:gen][g]["qmin"] <= qg[g in 1:G] <= ref[:gen][g]["qmax"])
    # Bi-directional branch flows
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= pf[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= qf[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    # 
    #   II. Constraints
    #
    # Nodal power balance
    JuMP.@constraint(model, 
        kirchhoff_active[i in 1:N],
        sum(pf[a] for a in ref[:bus_arcs][i]) ==
        sum(pg[g] for g in ref[:bus_gens][i]) -
        sum(load["pd"] for load in bus_loads[i]) -
        sum(shunt["gs"] for shunt in bus_shunts[i])*w[i]
    )
    JuMP.@constraint(model,
        kirchhoff_reactive[i in 1:N],
        sum(qf[a] for a in ref[:bus_arcs][i]) ==
        sum(qg[g] for g in ref[:bus_gens][i]) -
        sum(load["qd"] for load in bus_loads[i]) +
        sum(shunt["bs"] for shunt in bus_shunts[i])*w[i]
    )

    # Branch power flow physics and limit constraints
    # We pre-allocate constraint containers for simplicity
    model[:thermal_limit_fr] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:thermal_limit_to] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:voltage_difference_limit_ub] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:voltage_difference_limit_lb] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:ohm_active_fr] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:ohm_active_to] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:ohm_reactive_fr] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:ohm_reactive_to] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:jabr] = Vector{JuMP.ConstraintRef}(undef, E)
    for (i,branch) in ref[:branch]
        data["branch"]["$i"]["br_status"] == 0 && continue  # skip branch

        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        bp_idx = (branch["f_bus"], branch["t_bus"])

        p_fr = pf[f_idx]
        q_fr = qf[f_idx]
        p_to = pf[t_idx]
        q_to = qf[t_idx]

        w_fr = w[branch["f_bus"]]
        w_to = w[branch["t_bus"]]
        wr_br = wr[i]
        wi_br = wi[i]

        g, b = PM.calc_branch_y(branch)
        tr, ti = PM.calc_branch_t(branch)
        ttm = tr^2 + ti^2
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]

        # From side of the branch flow
        model[:ohm_active_fr][i] = JuMP.@constraint(model, p_fr ==  (g+g_fr)/ttm*w_fr + (-g*tr+b*ti)/ttm*(wr_br) + (-b*tr-g*ti)/ttm*(wi_br) )
        model[:ohm_reactive_fr][i] = JuMP.@constraint(model, q_fr == -(b+b_fr)/ttm*w_fr - (-b*tr-g*ti)/ttm*(wr_br) + (-g*tr+b*ti)/ttm*(wi_br) )

        # To side of the branch flow
        model[:ohm_active_to][i] = JuMP.@constraint(model, p_to ==  (g+g_to)*w_to + (-g*tr-b*ti)/ttm*(wr_br) + (-b*tr+g*ti)/ttm*(-wi_br) )
        model[:ohm_reactive_to][i] = JuMP.@constraint(model, q_to == -(b+b_to)*w_to - (-b*tr+g*ti)/ttm*(wr_br) + (-g*tr-b*ti)/ttm*(-wi_br) )

        # Voltage angle difference limit
        model[:voltage_difference_limit_ub][i] = JuMP.@constraint(model, wi_br <= tan(branch["angmax"])*wr_br)
        model[:voltage_difference_limit_lb][i] = JuMP.@constraint(model, wi_br >= tan(branch["angmin"])*wr_br)

        # Apparent power limit, from side and to side
        if OPF == PM.SOCWRPowerModel
            model[:thermal_limit_fr][i] = JuMP.@constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
            model[:thermal_limit_to][i] = JuMP.@constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
        elseif OPF == PM.SOCWRConicPowerModel
            model[:thermal_limit_fr][i] = JuMP.@constraint(model, [branch["rate_a"], p_fr, q_fr] in JuMP.SecondOrderCone())
            model[:thermal_limit_to][i] = JuMP.@constraint(model, [branch["rate_a"], p_to, q_to] in JuMP.SecondOrderCone())
        end

        # Jabr inequality (one per branch)
        if OPF == PM.SOCWRPowerModel
            model[:jabr][i] = JuMP.@constraint(model,
                wr[i]^2 + wi[i]^2 <= w[branch["f_bus"]]*w[branch["t_bus"]]
            )
        elseif OPF == PM.SOCWRConicPowerModel
            model[:jabr][i] = JuMP.@constraint(model,
                [w[branch["f_bus"]] / sqrt(2), w[branch["t_bus"]] / sqrt(2), wr[i], wi[i]] in JuMP.RotatedSecondOrderCone()
            )
        end
    end

    #
    #   III. Objective
    #
    l, u = extrema(gen["cost"][1] for (i, gen) in ref[:gen])
    (l == u == 0.0) || @warn "Data $(data["name"]) has quadratic cost terms; those terms are being ignored"
    JuMP.@objective(model, Min, sum(
        gen["cost"][2]*pg[i] + gen["cost"][3]
        for (i,gen) in ref[:gen]
    ))

    return OPFModel{OPF}(data, model)
end

"""
    _extract_solution(model, data)

Extract SOC-OPF solution from optimization model.
The model must have been solved before.
"""
function extract_result(opf::OPFModel{OPF}) where {OPF <: Union{PM.SOCWRPowerModel,PM.SOCWRConicPowerModel}}
    data  = opf.data
    model = opf.model

    # Pre-process data
    ref = PM.build_ref(data)[:it][:pm][:nw][0]
    N = length(ref[:bus])
    G = length(ref[:gen])
    E = length(data["branch"])

    # Build the solution dictionary
    res = Dict{String,Any}()
    res["opf_model"] = string(model.ext[:opf_model])
    res["objective"] = JuMP.objective_value(model)
    res["objective_lb"] = -Inf
    res["optimizer"] = JuMP.solver_name(model)
    res["solve_time"] = JuMP.solve_time(model)
    res["termination_status"] = JuMP.termination_status(model)
    res["primal_status"] = JuMP.primal_status(model)
    res["dual_status"] = JuMP.dual_status(model)
    res["solution"] = sol = Dict{String,Any}()

    sol["per_unit"] = get(data, "per_unit", false)
    sol["baseMVA"]  = get(data, "baseMVA", 100.0)

    ### populate branches, gens, buses ###

    sol["bus"] = Dict{String,Any}()
    sol["branch"] = Dict{String,Any}()
    sol["gen"] = Dict{String,Any}()

    for bus in 1:N
        sol["bus"]["$bus"] = Dict(
            "w" => value(model[:w][bus]),
            # dual vars
            "lam_kirchhoff_active" => dual(model[:kirchhoff_active][bus]),
            "lam_kirchhoff_reactive" => dual(model[:kirchhoff_reactive][bus]),
            "mu_w_lb" => dual(LowerBoundRef(model[:w][bus])),
            "mu_w_ub" => dual(UpperBoundRef(model[:w][bus])),
        )
    end

    wbounds = model.ext[:wbounds]
    for b in 1:E
        if data["branch"]["$b"]["br_status"] == 0
            # branch is under outage --> we set everything to zero
            # The variables `wr` and `wi` don't really have a meaning anymore,
            #   so we set them to zero by convention.
            sol["branch"]["$b"] = brsol = Dict(
                "pf" => 0.0,
                "pt" => 0.0,
                "qf" => 0.0,
                "qt" => 0.0,
                "wr" => 0.0,
                "wi" => 0.0,
                # dual vars
                "mu_va_diff_ub" => 0.0,
                "mu_va_diff_lb" => 0.0,
                "mu_wr_lb" => 0.0,
                "mu_wr_ub" => 0.0,
                "mu_wi_lb" => 0.0,
                "mu_wi_ub" => 0.0,
                "lam_ohm_active_fr" => 0.0,
                "lam_ohm_active_to" => 0.0,
                "lam_ohm_reactive_fr" => 0.0,
                "lam_ohm_reactive_to" => 0.0,
            )
            if OPF == PM.SOCWRPowerModel
                brsol["mu_sm_to"] = 0.0
                brsol["mu_sm_fr"] = 0.0
                brsol["mu_voltage_prod_quad"] = 0.0
            elseif OPF == PM.SOCWRConicPowerModel
                # conic duals (unrolled)
                brsol["nu_voltage_prod_soc_1"] = 0.0
                brsol["nu_voltage_prod_soc_2"] = 0.0
                brsol["nu_voltage_prod_soc_3"] = 0.0
                brsol["nu_voltage_prod_soc_4"] = 0.0
                brsol["nu_sm_to_1"] = 0.0
                brsol["nu_sm_to_2"] = 0.0
                brsol["nu_sm_to_3"] = 0.0
                brsol["nu_sm_fr_1"] = 0.0
                brsol["nu_sm_fr_2"] = 0.0
                brsol["nu_sm_fr_3"] = 0.0
            end
        else
            bp = (ref[:branch][b]["f_bus"], ref[:branch][b]["t_bus"])
            sol["branch"]["$b"] = brsol = Dict(
                "pf" => value(model[:pf][(b,ref[:branch][b]["f_bus"],ref[:branch][b]["t_bus"])]),
                "pt" => value(model[:pf][(b,ref[:branch][b]["t_bus"],ref[:branch][b]["f_bus"])]),
                "qf" => value(model[:qf][(b,ref[:branch][b]["f_bus"],ref[:branch][b]["t_bus"])]),
                "qt" => value(model[:qf][(b,ref[:branch][b]["t_bus"],ref[:branch][b]["f_bus"])]),
                "wr" => value(model[:wr][b]),
                "wi" => value(model[:wi][b]),
                # dual vars
                "mu_va_diff_ub" => dual(model[:voltage_difference_limit_ub][b]),
                "mu_va_diff_lb" => dual(model[:voltage_difference_limit_lb][b]),
                "mu_wr_lb" => wbounds ? dual(LowerBoundRef(model[:wr][b])) : 0.0,
                "mu_wr_ub" => wbounds ? dual(UpperBoundRef(model[:wr][b])) : 0.0,
                "mu_wi_lb" => wbounds ? dual(LowerBoundRef(model[:wi][b])) : 0.0,
                "mu_wi_ub" => wbounds ? dual(UpperBoundRef(model[:wi][b])) : 0.0,
                "lam_ohm_active_fr" => dual(model[:ohm_active_fr][b]),
                "lam_ohm_active_to" => dual(model[:ohm_active_to][b]),
                "lam_ohm_reactive_fr" => dual(model[:ohm_reactive_fr][b]),
                "lam_ohm_reactive_to" => dual(model[:ohm_reactive_to][b]),
            )
            
            # The following dual variables depend on whether we have a QCP or SOC form,
            #   so they are handled separately
            if OPF == PM.SOCWRPowerModel
                brsol["mu_sm_to"] = dual(model[:thermal_limit_to][b])
                brsol["mu_sm_fr"] = dual(model[:thermal_limit_fr][b])
                brsol["mu_voltage_prod_quad"] = dual(model[:jabr][b])
            elseif OPF == PM.SOCWRConicPowerModel
                nu_voltage_prod_soc = dual(model[:jabr][b])
                brsol["nu_voltage_prod_soc_1"] = nu_voltage_prod_soc[1]
                brsol["nu_voltage_prod_soc_2"] = nu_voltage_prod_soc[2]
                brsol["nu_voltage_prod_soc_3"] = nu_voltage_prod_soc[3]
                brsol["nu_voltage_prod_soc_4"] = nu_voltage_prod_soc[4]
                nu_sm_to = dual(model[:thermal_limit_to][b])
                brsol["nu_sm_to_1"] = nu_sm_to[1]
                brsol["nu_sm_to_2"] = nu_sm_to[2]
                brsol["nu_sm_to_3"] = nu_sm_to[3]
                nu_sm_fr = dual(model[:thermal_limit_fr][b])
                brsol["nu_sm_fr_1"] = nu_sm_fr[1]
                brsol["nu_sm_fr_2"] = nu_sm_fr[2]
                brsol["nu_sm_fr_3"] = nu_sm_fr[3]
            end
        end
    end
    
    for g in 1:G
        sol["gen"]["$g"] = Dict(
            "pg" => value(model[:pg][g]),
            "qg" => value(model[:qg][g]),
            # dual vars
            "mu_pg_lb" => dual(LowerBoundRef(model[:pg][g])),
            "mu_pg_ub" => dual(UpperBoundRef(model[:pg][g])),
            "mu_qg_lb" => dual(LowerBoundRef(model[:qg][g])),
            "mu_qg_ub" => dual(UpperBoundRef(model[:qg][g])),
        )
    end
    
    return res
end

function json2h5(::Type{OPF}, res) where{OPF <: Union{PM.SOCWRPowerModel,PM.SOCWRConicPowerModel}}
    sol = res["solution"]
    N = length(sol["bus"])
    E = length(sol["branch"])
    G = length(sol["gen"])

    res_h5 = Dict{String,Any}(
        "meta" => Dict{String,Any}(
            "termination_status" => res["termination_status"],
            "primal_status" => res["primal_status"],
            "dual_status" => res["dual_status"],
            "solve_time" => res["solve_time"],
        ),
    )

    res_h5["primal"] = pres_h5 = Dict{String,Any}(
        "w"  => zeros(Float64, N),
        "pg" => zeros(Float64, G),
        "qg" => zeros(Float64, G),
        "pf" => zeros(Float64, E),
        "qf" => zeros(Float64, E),
        "pt" => zeros(Float64, E),
        "qt" => zeros(Float64, E),
        "wr" => zeros(Float64, E),
        "wi" => zeros(Float64, E),
    )
    res_h5["dual"] = dres_h5 = Dict{String,Any}(
        "mu_w_lb"                => zeros(Float64, N),
        "mu_w_ub"                => zeros(Float64, N),
        "lam_kirchhoff_active"   => zeros(Float64, N),
        "lam_kirchhoff_reactive" => zeros(Float64, N),
        "mu_pg_lb"               => zeros(Float64, G),
        "mu_pg_ub"               => zeros(Float64, G),
        "mu_qg_lb"               => zeros(Float64, G),
        "mu_qg_ub"               => zeros(Float64, G),
        # branhc dual vars
        "mu_va_diff_lb"          => zeros(Float64, E),
        "mu_va_diff_ub"          => zeros(Float64, E),
        "mu_wr_lb"               => zeros(Float64, E),
        "mu_wr_ub"               => zeros(Float64, E),
        "mu_wi_lb"               => zeros(Float64, E),
        "mu_wi_ub"               => zeros(Float64, E),
        "lam_ohm_active_fr"      => zeros(Float64, E),
        "lam_ohm_active_to"      => zeros(Float64, E),
        "lam_ohm_reactive_fr"    => zeros(Float64, E),
        "lam_ohm_reactive_to"    => zeros(Float64, E),
    )
    # The following dual variables depend on whether we have a SOC or QCP form
    if OPF == PM.SOCWRPowerModel
        dres_h5["mu_sm_to"] = zeros(Float64, E)
        dres_h5["mu_sm_fr"] = zeros(Float64, E)
        dres_h5["mu_voltage_prod_quad"] = zeros(Float64, E)
    elseif OPF == PM.SOCWRConicPowerModel
        # conic duals (unrolled)
        dres_h5["nu_voltage_prod_soc_1"] = zeros(Float64, E)
        dres_h5["nu_voltage_prod_soc_2"] = zeros(Float64, E)
        dres_h5["nu_voltage_prod_soc_3"] = zeros(Float64, E)
        dres_h5["nu_voltage_prod_soc_4"] = zeros(Float64, E)
        dres_h5["nu_sm_to_1"] = zeros(Float64, E)
        dres_h5["nu_sm_to_2"] = zeros(Float64, E)
        dres_h5["nu_sm_to_3"] = zeros(Float64, E)
        dres_h5["nu_sm_fr_1"] = zeros(Float64, E)
        dres_h5["nu_sm_fr_2"] = zeros(Float64, E)
        dres_h5["nu_sm_fr_3"] = zeros(Float64, E)
    end

    # extract from ACOPF solution
    for i in 1:N
        bsol = sol["bus"]["$i"]

        pres_h5["w"][i] = bsol["w"]

        dres_h5["mu_w_lb"][i] = bsol["mu_w_lb"]
        dres_h5["mu_w_ub"][i] = bsol["mu_w_ub"]
        dres_h5["lam_kirchhoff_active"][i] = bsol["lam_kirchhoff_active"]
        dres_h5["lam_kirchhoff_reactive"][i] = bsol["lam_kirchhoff_reactive"]
    end
    for g in 1:G
        gsol = sol["gen"]["$g"]

        pres_h5["pg"][g] = gsol["pg"]
        pres_h5["qg"][g] = gsol["qg"]

        dres_h5["mu_pg_lb"][g] = gsol["mu_pg_lb"]
        dres_h5["mu_pg_ub"][g] = gsol["mu_pg_ub"]
        dres_h5["mu_qg_lb"][g] = gsol["mu_qg_lb"]
        dres_h5["mu_qg_ub"][g] = gsol["mu_qg_ub"]
    end
    for e in 1:E
        brsol = sol["branch"]["$e"]

        for pvar in ["pf", "qf", "pt", "qt", "wr", "wi"]
            pres_h5[pvar][e] = brsol[pvar]
        end
        for dvar in [
            "mu_va_diff_lb",
            "mu_va_diff_ub",
            "mu_wr_lb",
            "mu_wr_ub",
            "mu_wi_lb",
            "mu_wi_ub",
            "lam_ohm_active_fr",
            "lam_ohm_active_to",
            "lam_ohm_reactive_fr",
            "lam_ohm_reactive_to",
        ]
            dres_h5[dvar][e] = brsol[dvar]
        end
        
        # The following dual variables depend on whether we have a SOC or QCP form
        if OPF == PM.SOCWRPowerModel
            dres_h5["mu_sm_to"][e] = brsol["mu_sm_to"]
            dres_h5["mu_sm_fr"][e] = brsol["mu_sm_fr"]
            dres_h5["mu_voltage_prod_quad"][e] = brsol["mu_voltage_prod_quad"]
        elseif OPF == PM.SOCWRConicPowerModel
            # conic duals (unrolled)
            dres_h5["nu_voltage_prod_soc_1"][e] = brsol["nu_voltage_prod_soc_1"]
            dres_h5["nu_voltage_prod_soc_2"][e] = brsol["nu_voltage_prod_soc_2"]
            dres_h5["nu_voltage_prod_soc_3"][e] = brsol["nu_voltage_prod_soc_3"]
            dres_h5["nu_voltage_prod_soc_4"][e] = brsol["nu_voltage_prod_soc_4"]
            dres_h5["nu_sm_to_1"][e] = brsol["nu_sm_to_1"]
            dres_h5["nu_sm_to_2"][e] = brsol["nu_sm_to_2"]
            dres_h5["nu_sm_to_3"][e] = brsol["nu_sm_to_3"]
            dres_h5["nu_sm_fr_1"][e] = brsol["nu_sm_fr_1"]
            dres_h5["nu_sm_fr_2"][e] = brsol["nu_sm_fr_2"]
            dres_h5["nu_sm_fr_3"][e] = brsol["nu_sm_fr_3"]
        end
    end

    return res_h5
end
