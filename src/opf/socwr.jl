"""
    build_soc_opf(data, optimizer)

Build an SOC-OPF model.

This implementation is based on the SOC-OPF formulation of PMAnnex.jl
    https://github.com/lanl-ansi/PMAnnex.jl/blob/f303f3c3c61e2d1a050ee7651fa6e8abc4055b55/src/model/opf.jl
"""
function build_opf(::Type{PM.SOCWRPowerModel}, data::Dict{String,Any}, optimizer)
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

    model = JuMP.Model(optimizer)
    model.ext[:opf_model] = PM.SOCWRPowerModel

    #
    #   I. Variables
    #
    # nodal voltage
    JuMP.@variable(model, ref[:bus][i]["vmin"]^2 <= w[i in 1:N] <= ref[:bus][i]["vmax"]^2, start=1.001)

    wr_min, wr_max, wi_min, wi_max = PM.ref_calc_voltage_product_bounds(ref[:buspairs])

    JuMP.@variable(model, wr_min[bp] <= wr[bp in keys(ref[:buspairs])] <= wr_max[bp], start=1.0)
    JuMP.@variable(model, wi_min[bp] <= wi[bp in keys(ref[:buspairs])] <= wi_max[bp])
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
        wr_br = wr[bp_idx]
        wi_br = wi[bp_idx]

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
        model[:thermal_limit_fr][i] = JuMP.@constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        model[:thermal_limit_to][i] = JuMP.@constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
    end
    
    # Voltage product relaxation (quadratic form)
    JuMP.@constraint(model,
        voltage_prod_quadratic[(i, j) in keys(ref[:buspairs])],
        wr[(i,j)]^2 + wi[(i,j)]^2 <= w[i]*w[j]
    )

    #
    #   III. Objective
    #
    JuMP.@objective(model, Min, sum(
        gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3]
        for (i,gen) in ref[:gen]
    ))

    return OPFModel{PM.SOCWRPowerModel}(data, model)
end

"""
    _extract_solution(model, data)

Extract SOC-OPF solution from optimization model.
The model must have been solved before.
"""
function extract_result(opf::OPFModel{PM.SOCWRPowerModel})
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
            "lam_pb_active" => dual(model[:kirchhoff_active][bus]),
            "lam_pb_reactive" => dual(model[:kirchhoff_reactive][bus]),
            "mu_w_lb" => dual(LowerBoundRef(model[:w][bus])),
            "mu_w_ub" => dual(UpperBoundRef(model[:w][bus])),
        )
    end

    for b in 1:E
        if data["branch"]["$b"]["br_status"] == 0
            # branch is under outage --> we set everything to zero
            # The variables `wr` and `wi` don't really have a meaning anymore,
            #   so we set them to zero by convention.
            sol["branch"]["$b"] = Dict(
                "pf" => 0.0,
                "pt" => 0.0,
                "qf" => 0.0,
                "qt" => 0.0,
                "wr" => 0.0,
                "wi" => 0.0,
                # dual vars
                "mu_sm_to" => 0.0,
                "mu_sm_fr" => 0.0,
                "mu_va_diff_ub" => 0.0,
                "mu_va_diff_lb" => 0.0,
                "mu_voltage_prod_quad" => 0.0,
                "mu_wr_lb" => 0.0,
                "mu_wr_ub" => 0.0,
                "mu_wi_lb" => 0.0,
                "mu_wi_ub" => 0.0,
                "lam_ohm_active_fr" => 0.0,
                "lam_ohm_active_to" => 0.0,
                "lam_ohm_reactive_fr" => 0.0,
                "lam_ohm_reactive_to" => 0.0,
            )
        else
            bp = (ref[:branch][b]["f_bus"], ref[:branch][b]["t_bus"])
            sol["branch"]["$b"] = Dict(
                "pf" => value(model[:pf][(b,ref[:branch][b]["f_bus"],ref[:branch][b]["t_bus"])]),
                "pt" => value(model[:pf][(b,ref[:branch][b]["t_bus"],ref[:branch][b]["f_bus"])]),
                "qf" => value(model[:qf][(b,ref[:branch][b]["f_bus"],ref[:branch][b]["t_bus"])]),
                "qt" => value(model[:qf][(b,ref[:branch][b]["t_bus"],ref[:branch][b]["f_bus"])]),
                "wr" => value(model[:wr][bp]),
                "wi" => value(model[:wi][bp]),
                # dual vars
                "mu_sm_to" => dual(model[:thermal_limit_to][b]),
                "mu_sm_fr" => dual(model[:thermal_limit_fr][b]),
                "mu_va_diff_ub" => dual(model[:voltage_difference_limit_ub][b]),
                "mu_va_diff_lb" => dual(model[:voltage_difference_limit_lb][b]),
                "mu_voltage_prod_quad" => dual(model[:voltage_prod_quadratic][bp]),
                "mu_wr_lb" => dual(LowerBoundRef(model[:wr][bp])),
                "mu_wr_ub" => dual(UpperBoundRef(model[:wr][bp])),
                "mu_wi_lb" => dual(LowerBoundRef(model[:wi][bp])),
                "mu_wi_ub" => dual(UpperBoundRef(model[:wi][bp])),
                "lam_ohm_active_fr" => dual(model[:ohm_active_fr][b]),
                "lam_ohm_active_to" => dual(model[:ohm_active_to][b]),
                "lam_ohm_reactive_fr" => dual(model[:ohm_reactive_fr][b]),
                "lam_ohm_reactive_to" => dual(model[:ohm_reactive_to][b]),
            )
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
