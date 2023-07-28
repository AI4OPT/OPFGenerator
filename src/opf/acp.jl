"""
    build_acopf(data, optimizer)

Build an AC-OPF model.

This implementation is based on the AC-OPF formulation of Rosetta-OPF
    https://github.com/lanl-ansi/rosetta-opf/blob/38a951326df3156d79dcdc49c8010aa29905b05d/jump.jl
"""
function build_opf(::Type{PM.ACPPowerModel}, data::Dict{String,Any}, optimizer)
    # Cleanup and pre-process data
    PM.standardize_cost_terms!(data, order=2)
    PM.calc_thermal_limits!(data)
    ref = PM.build_ref(data)[:it][:pm][:nw][0]

    # Check that slack bus is unique
    length(ref[:ref_buses]) == 1 || error("Data has $(length(ref[:ref_buses])) slack buses (expected 1).")
    i0 = first(keys(ref[:ref_buses]))
    
    # Grab some data
    N = length(ref[:bus])
    G = length(ref[:gen])
    E = length(data["branch"])
    L = length(ref[:load])
    bus_loads = [
        [ref[:load][l] for l in ref[:bus_loads][i]]
        for i in 1:N
    ]
    bus_shunts = [
        [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        for i in 1:N
    ]

    model = JuMP.Model(optimizer)
    model.ext[:opf_model] = PM.ACPPowerModel  # for internal checks

    #
    #   I. Variables
    #
    # nodal voltage
    JuMP.@variable(model, va[1:N])
    JuMP.@variable(model, ref[:bus][i]["vmin"] <= vm[i in 1:N] <= ref[:bus][i]["vmax"], start=1.0)
    # Active and reactive dispatch
    JuMP.@variable(model, ref[:gen][g]["pmin"] <= pg[g in 1:G] <= ref[:gen][g]["pmax"])
    JuMP.@variable(model, ref[:gen][g]["qmin"] <= qg[g in 1:G] <= ref[:gen][g]["qmax"])
    # Bi-directional branch flows
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= pf[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= qf[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    # 
    #   II. Constraints
    #
    # Slack bus
    JuMP.@constraint(model, slack_bus, va[i0] == 0.0)

    # Nodal power balance
    JuMP.@constraint(model, 
        kirchhoff_active[i in 1:N],
        sum(pf[a] for a in ref[:bus_arcs][i]) ==
        sum(pg[g] for g in ref[:bus_gens][i]) -
        sum(load["pd"] for load in bus_loads[i]) -
        sum(shunt["gs"] for shunt in bus_shunts[i])*vm[i]^2
    )
    JuMP.@constraint(model,
        kirchhoff_reactive[i in 1:N],
        sum(qf[a] for a in ref[:bus_arcs][i]) ==
        sum(qg[g] for g in ref[:bus_gens][i]) -
        sum(load["qd"] for load in bus_loads[i]) +
        sum(shunt["bs"] for shunt in bus_shunts[i])*vm[i]^2
    )

    # Branch power flow physics and limit constraints
    # We pre-allocate constraint containers for simplicity
    model[:thermal_limit_fr] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:thermal_limit_to] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:voltage_difference_limit] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:ohm_active_fr] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:ohm_active_to] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:ohm_reactive_fr] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:ohm_reactive_to] = Vector{JuMP.ConstraintRef}(undef, E)
    for (i,branch) in ref[:branch]
        data["branch"]["$i"]["br_status"] == 0 && continue  # skip branch

        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = pf[f_idx]
        q_fr = qf[f_idx]
        p_to = pf[t_idx]
        q_to = qf[t_idx]

        vm_fr = vm[branch["f_bus"]]
        vm_to = vm[branch["t_bus"]]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PM.calc_branch_y(branch)
        tr, ti = PM.calc_branch_t(branch)
        ttm = tr^2 + ti^2
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]

        # From side of the branch flow
        model[:ohm_active_fr][i] = JuMP.@NLconstraint(model, p_fr ==  (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        model[:ohm_reactive_fr][i] = JuMP.@NLconstraint(model, q_fr == -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        model[:ohm_active_to][i] = JuMP.@NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        model[:ohm_reactive_to][i] = JuMP.@NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Voltage angle difference limit
        model[:voltage_difference_limit][i] = JuMP.@constraint(model, branch["angmin"] <= va_fr - va_to <= branch["angmax"])

        # Apparent power limit, from side and to side
        model[:thermal_limit_fr][i] = JuMP.@constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        model[:thermal_limit_to][i] = JuMP.@constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
    end

    #
    #   III. Objective
    #
    # check that we don't have quadratic objective coeffs
    l, u = extrema(gen["cost"][1] for (i, gen) in ref[:gen])
    (l == u == 0.0) || @warn "Data $(data["name"]) has quadratic cost terms; those terms are being ignored"
    JuMP.@objective(model, Min, sum(
        gen["cost"][2]*pg[i] + gen["cost"][3]
        for (i,gen) in ref[:gen]
    ))

    return OPFModel{PM.ACPPowerModel}(data, model)
end

"""
    _extract_acopf_solution(model, data)

Extract ACOPF solution from optimization model.
The model must have been solved before.
"""
function extract_result(opf::OPFModel{PM.ACPPowerModel})
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
            "vm" => value(model[:vm][bus]),
            "va" => value(model[:va][bus]),
            # dual vars
            "lam_pb_active" => dual(model[:kirchhoff_active][bus]),
            "lam_pb_reactive" => dual(model[:kirchhoff_reactive][bus]),
            "mu_vm_lb" => dual(LowerBoundRef(model[:vm][bus])),
            "mu_vm_ub" => dual(UpperBoundRef(model[:vm][bus]))
        )
    end

    for b in 1:E
        if data["branch"]["$b"]["br_status"] == 0
            # branch is under outage --> we set everything to zero
            sol["branch"]["$b"] = Dict(
                "pf" => 0.0,
                "pt" => 0.0,
                "qf" => 0.0,
                "qt" => 0.0,
                # dual vars
                "mu_sm_to" => 0.0,
                "mu_sm_fr" => 0.0,
                "mu_va_diff" => 0.0,
                "lam_ohm_active_fr" => 0.0,
                "lam_ohm_active_to" => 0.0,
                "lam_ohm_reactive_fr" => 0.0,
                "lam_ohm_reactive_to" => 0.0,
            )
        else
            sol["branch"]["$b"] = Dict(
                "pf" => value(model[:pf][(b,ref[:branch][b]["f_bus"],ref[:branch][b]["t_bus"])]),
                "pt" => value(model[:pf][(b,ref[:branch][b]["t_bus"],ref[:branch][b]["f_bus"])]),
                "qf" => value(model[:qf][(b,ref[:branch][b]["f_bus"],ref[:branch][b]["t_bus"])]),
                "qt" => value(model[:qf][(b,ref[:branch][b]["t_bus"],ref[:branch][b]["f_bus"])]),
                # dual vars
                "mu_sm_to" => dual(model[:thermal_limit_to][b]),
                "mu_sm_fr" => dual(model[:thermal_limit_fr][b]),
                "mu_va_diff" => dual(model[:voltage_difference_limit][b]),
                "lam_ohm_active_fr" => dual(model[:ohm_active_fr][b]),
                "lam_ohm_active_to" => dual(model[:ohm_active_to][b]),
                "lam_ohm_reactive_fr" => dual(model[:ohm_reactive_fr][b]),
                "lam_ohm_reactive_to" => dual(model[:ohm_reactive_to][b])
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
            "mu_qg_ub" => dual(UpperBoundRef(model[:qg][g]))
        )
    end 
    
    return res
end
