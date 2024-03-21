"""
    build_acopf(data, optimizer)

Build an AC-OPF model.

This implementation is based on the AC-OPF formulation of Rosetta-OPF
    https://github.com/lanl-ansi/rosetta-opf/blob/38a951326df3156d79dcdc49c8010aa29905b05d/jump.jl
"""
function build_opf(::Type{PM.ACPPowerModel}, data::Dict{String,Any}, optimizer;
    T=Float64,    
)
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

    model = JuMP.GenericModel{T}(optimizer)
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
        model[:ohm_active_fr][i] = JuMP.@constraint(model, p_fr ==  (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        model[:ohm_reactive_fr][i] = JuMP.@constraint(model, q_fr == -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        model[:ohm_active_to][i] = JuMP.@constraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        model[:ohm_reactive_to][i] = JuMP.@constraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )

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

    sol["singleton"] = Dict(
        "lam_slack_bus" => dual(model[:slack_bus]),
    )

    return res
end

function json2h5(::Type{PM.ACPPowerModel}, res)
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
        "vm" => zeros(Float64, N),
        "va" => zeros(Float64, N),
        "pg" => zeros(Float64, G),
        "qg" => zeros(Float64, G),
        "pf" => zeros(Float64, E),
        "qf" => zeros(Float64, E),
        "pt" => zeros(Float64, E),
        "qt" => zeros(Float64, E),
    )
    res_h5["dual"] = dres_h5 = Dict{String,Any}(
        "mu_vm_lb"               => zeros(Float64, N),
        "mu_vm_ub"               => zeros(Float64, N),
        "lam_kirchhoff_active"   => zeros(Float64, N),
        "lam_kirchhoff_reactive" => zeros(Float64, N),
        "mu_pg_lb"               => zeros(Float64, G),
        "mu_pg_ub"               => zeros(Float64, G),
        "mu_qg_lb"               => zeros(Float64, G),
        "mu_qg_ub"               => zeros(Float64, G),
        "mu_sm_fr"               => zeros(Float64, E),
        "mu_sm_to"               => zeros(Float64, E),
        "lam_ohm_active_fr"      => zeros(Float64, E),
        "lam_ohm_active_to"      => zeros(Float64, E),
        "lam_ohm_reactive_fr"    => zeros(Float64, E),
        "lam_ohm_reactive_to"    => zeros(Float64, E),
        "mu_va_diff"             => zeros(Float64, E),
        "lam_slack_bus"          => 0.0,
    )

    # extract from ACOPF solution
    for i in 1:N
        bsol = sol["bus"]["$i"]

        pres_h5["vm"][i] = bsol["vm"]
        pres_h5["va"][i] = bsol["va"]

        dres_h5["mu_vm_lb"][i] = bsol["mu_vm_lb"]
        dres_h5["mu_vm_ub"][i] = bsol["mu_vm_ub"]
        dres_h5["lam_kirchhoff_active"][i] = bsol["lam_pb_active"]
        dres_h5["lam_kirchhoff_reactive"][i] = bsol["lam_pb_reactive"]
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

        pres_h5["pf"][e] = brsol["pf"]
        pres_h5["qf"][e] = brsol["qf"]
        pres_h5["pt"][e] = brsol["pt"]
        pres_h5["qt"][e] = brsol["qt"]
        
        dres_h5["mu_sm_fr"][e] = brsol["mu_sm_fr"]
        dres_h5["mu_sm_to"][e] = brsol["mu_sm_to"]
        dres_h5["lam_ohm_active_fr"][e] = brsol["lam_ohm_active_fr"]
        dres_h5["lam_ohm_active_to"][e] = brsol["lam_ohm_active_to"]
        dres_h5["lam_ohm_reactive_fr"][e] = brsol["lam_ohm_reactive_fr"]
        dres_h5["lam_ohm_reactive_to"][e] = brsol["lam_ohm_reactive_to"]
        dres_h5["mu_va_diff"][e] = brsol["mu_va_diff"]
    end

    dres_h5["lam_slack_bus"] = sol["singleton"]["lam_slack_bus"]

    return res_h5
end
