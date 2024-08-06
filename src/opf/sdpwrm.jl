"""
    build_opf(data, optimizer)

Build an SDP-OPF model.
"""
function build_opf(::Type{PM.SDPWRMPowerModel}, data::Dict{String,Any}, optimizer;
    T=Float64,
)
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

    model = JuMP.GenericModel{T}(optimizer)
    model.ext[:opf_model] = PM.SDPWRMPowerModel

    #
    #   I. Variables
    #
    # nodal voltage
    WR_start = zeros(N, N) + I
    JuMP.@variable(model, WR[i=1:N, j=1:N], Symmetric, start=WR_start[i,j])
    JuMP.@variable(model, WI[1:N, 1:N] in SkewSymmetricMatrixSpace(), start=0.0)

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
        sum(pg[g] for g in ref[:bus_gens][i])
        - sum(pf[a] for a in ref[:bus_arcs][i])
        - sum(shunt["gs"] for shunt in bus_shunts[i]) * WR[i, i]
        == sum(load["pd"] for load in bus_loads[i])
    )
    JuMP.@constraint(model,
        kirchhoff_reactive[i in 1:N],
        sum(qg[g] for g in ref[:bus_gens][i]) 
        - sum(qf[a] for a in ref[:bus_arcs][i])
        + sum(shunt["bs"] for shunt in bus_shunts[i]) * WR[i, i]
        ==
        sum(load["qd"] for load in bus_loads[i])
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

        w_fr = WR[branch["f_bus"], branch["f_bus"]]
        w_to = WR[branch["t_bus"], branch["t_bus"]]
        wr_br = WR[branch["f_bus"], branch["t_bus"]]
        wi_br = WI[branch["f_bus"], branch["t_bus"]]

        g, b = PM.calc_branch_y(branch)
        tr, ti = PM.calc_branch_t(branch)
        ttm = tr^2 + ti^2
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]

        # From side of the branch flow
        model[:ohm_active_fr][i] = JuMP.@constraint(model,
              ((g+g_fr)/ttm) * w_fr 
            + ((-g*tr+b*ti)/ttm) * wr_br 
            + ((-b*tr-g*ti)/ttm) * wi_br
            - p_fr
            ==
            0
        )
        model[:ohm_reactive_fr][i] = JuMP.@constraint(model,
            - ((b+b_fr)/ttm) * w_fr 
            - ((-b*tr-g*ti)/ttm) * wr_br
            + ((-g*tr+b*ti)/ttm) * wi_br
            - q_fr
            ==
            0
        )

        # To side of the branch flow
        model[:ohm_active_to][i] = JuMP.@constraint(model,
              ((g+g_to)) * w_to 
            + ((-g*tr-b*ti)/ttm) * wr_br 
            - ((-b*tr+g*ti)/ttm) * wi_br
            - p_to
            ==
            0
        )
        model[:ohm_reactive_to][i] = JuMP.@constraint(model,
            - ((b+b_to)) * w_to 
            - ((-b*tr+g*ti)/ttm) * wr_br
            - ((-g*tr-b*ti)/ttm) * wi_br 
            - q_to
            ==
            0
        )

        # Voltage angle difference limit
        model[:voltage_difference_limit_ub][i] = JuMP.@constraint(model, wi_br <= tan(branch["angmax"])*wr_br)
        model[:voltage_difference_limit_lb][i] = JuMP.@constraint(model, wi_br >= tan(branch["angmin"])*wr_br)

        # Apparent power limit, from side and to side
        model[:thermal_limit_fr][i] = JuMP.@constraint(model, [branch["rate_a"], p_fr, q_fr] in JuMP.SecondOrderCone())
        model[:thermal_limit_to][i] = JuMP.@constraint(model, [branch["rate_a"], p_to, q_to] in JuMP.SecondOrderCone())
    end

    # Nodal voltage bounds
    for (i, bus) in ref[:bus]
        wr_ii = WR[i, i]
        JuMP.set_lower_bound(wr_ii, (bus["vmin"])^2)
        JuMP.set_upper_bound(wr_ii, (bus["vmax"])^2)
    end

    # PSD constraint
    @constraint(model, S, [WR WI; -WI WR] in PSDCone())
    
    #
    #   III. Objective
    #
    l, u = extrema(gen["cost"][1] for (i, gen) in ref[:gen])
    (l == u == 0.0) || @warn "Data $(data["name"]) has quadratic cost terms; those terms are being ignored"
    JuMP.@objective(model, Min, sum(
        gen["cost"][2]*pg[i] + gen["cost"][3]
        for (i,gen) in ref[:gen]
    ))

    return OPFModel{PM.SDPWRMPowerModel}(data, model)
end

function update!(opf::OPFModel{PM.SDPWRMPowerModel}, data::Dict{String,Any})
    PM.standardize_cost_terms!(data, order=2)
    PM.calc_thermal_limits!(data)
    ref = PM.build_ref(data)[:it][:pm][:nw][0]

    opf.data = data

    N = length(ref[:bus])

    pd = [sum(ref[:load][l]["pd"] for l in ref[:bus_loads][i]; init=0.0) for i in 1:N]
    qd = [sum(ref[:load][l]["qd"] for l in ref[:bus_loads][i]; init=0.0) for i in 1:N]

    JuMP.set_normalized_rhs.(opf.model[:kirchhoff_active], pd)
    JuMP.set_normalized_rhs.(opf.model[:kirchhoff_reactive], qd)

    return nothing
end

"""
    _extract_result(model, data)

Extract SDP-OPF solution from optimization model.
The model must have been solved before.
"""
function extract_result(opf::OPFModel{PM.SDPWRMPowerModel})
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
    res["objective_lb"] = try JuMP.dual_objective_value(model) catch; NaN end
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
    sol["buspairs"] = Dict{String,Any}()

    for bus in 1:N
        sol["bus"]["$bus"] = Dict(
            "wm" => value(model[:WR][bus, bus]),
            # dual vars
            "lam_kirchhoff_active" => dual(model[:kirchhoff_active][bus]),
            "lam_kirchhoff_reactive" => dual(model[:kirchhoff_reactive][bus]),
            "mu_wm_lb" => dual(LowerBoundRef(model[:WR][bus, bus])),
            "mu_wm_ub" => dual(UpperBoundRef(model[:WR][bus, bus])),
        )
    end

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
                # dual vars
                "mu_va_diff_ub" => 0.0,
                "mu_va_diff_lb" => 0.0,
                "lam_ohm_active_fr" => 0.0,
                "lam_ohm_active_to" => 0.0,
                "lam_ohm_reactive_fr" => 0.0,
                "lam_ohm_reactive_to" => 0.0,
                # conic duals (vector-shaped)
                "nu_sm_to" => [0.0, 0.0, 0.0],
                "nu_sm_fr" => [0.0, 0.0, 0.0],
            )
        else
            br_f_bus = data["branch"]["$b"]["f_bus"]
            br_t_bus = data["branch"]["$b"]["t_bus"]
            sol["branch"]["$b"] = brsol = Dict{String,Any}(
                "pf" => value(model[:pf][(b,ref[:branch][b]["f_bus"],ref[:branch][b]["t_bus"])]),
                "pt" => value(model[:pf][(b,ref[:branch][b]["t_bus"],ref[:branch][b]["f_bus"])]),
                "qf" => value(model[:qf][(b,ref[:branch][b]["f_bus"],ref[:branch][b]["t_bus"])]),
                "qt" => value(model[:qf][(b,ref[:branch][b]["t_bus"],ref[:branch][b]["f_bus"])]),
                # dual vars
                "mu_va_diff_ub" => dual(model[:voltage_difference_limit_ub][b]),
                "mu_va_diff_lb" => dual(model[:voltage_difference_limit_lb][b]),
                "lam_ohm_active_fr" => dual(model[:ohm_active_fr][b]),
                "lam_ohm_active_to" => dual(model[:ohm_active_to][b]),
                "lam_ohm_reactive_fr" => dual(model[:ohm_reactive_fr][b]),
                "lam_ohm_reactive_to" => dual(model[:ohm_reactive_to][b]),
                "nu_sm_to" => dual(model[:thermal_limit_to][b]),
                "nu_sm_fr" => dual(model[:thermal_limit_fr][b]),
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

    for (p, _) in ref[:buspairs]
        sol["buspairs"]["$p"] = Dict(
            "wr" => value(model[:WR][p...]),
            "wi" => value(model[:WI][p...]),
        )
    end

    # TODO: save S sparsely
    sol["S"] = dual(model[:S])
    
    return res
end

function json2h5(::Type{PM.SDPWRMPowerModel}, res)
    sol = res["solution"]
    N = length(sol["bus"])
    E = length(sol["branch"])
    G = length(sol["gen"])

    res_h5 = Dict{String,Any}(
        "meta" => Dict{String,Any}(
            "termination_status" => string(res["termination_status"]),
            "primal_status" => string(res["primal_status"]),
            "dual_status" => string(res["dual_status"]),
            "solve_time" => res["solve_time"],
            "primal_objective_value" => res["objective"],
            "dual_objective_value" => res["objective_lb"],
        ),
    )

    res_h5["primal"] = pres_h5 = Dict{String,Any}(
        "wm" => zeros(Float64, N),
        "pg" => zeros(Float64, G),
        "qg" => zeros(Float64, G),
        "pf" => zeros(Float64, E),
        "qf" => zeros(Float64, E),
        "pt" => zeros(Float64, E),
        "qt" => zeros(Float64, E),
        "wr" => zeros(Float64, P),
        "wi" => zeros(Float64, P),
    )
    res_h5["dual"] = dres_h5 = Dict{String,Any}(
        "mu_wm_lb"               => zeros(Float64, N),
        "mu_wm_ub"               => zeros(Float64, N),
        "lam_kirchhoff_active"   => zeros(Float64, N),
        "lam_kirchhoff_reactive" => zeros(Float64, N),
        "mu_pg_lb"               => zeros(Float64, G),
        "mu_pg_ub"               => zeros(Float64, G),
        "mu_qg_lb"               => zeros(Float64, G),
        "mu_qg_ub"               => zeros(Float64, G),
        # branch dual vars
        "mu_va_diff_lb"          => zeros(Float64, E),
        "mu_va_diff_ub"          => zeros(Float64, E),
        "lam_ohm_active_fr"      => zeros(Float64, E),
        "lam_ohm_active_to"      => zeros(Float64, E),
        "lam_ohm_reactive_fr"    => zeros(Float64, E),
        "lam_ohm_reactive_to"    => zeros(Float64, E),
        # conic duals (unrolled)
        "nu_sm_to"               => zeros(Float64, E, 3),
        "nu_sm_fr"               => zeros(Float64, E, 3),
        "S"                      => zeros(Float64, 2*N, 2*N)
    )

    # extract from ACOPF solution
    for i in 1:N
        bsol = sol["bus"]["$i"]

        pres_h5["wm"][i] = bsol["wm"]

        dres_h5["mu_wm_lb"][i] = bsol["mu_wm_lb"]
        dres_h5["mu_wm_ub"][i] = bsol["mu_wm_ub"]
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

        for pvar in ["pf", "qf", "pt", "qt"]
            pres_h5[pvar][e] = brsol[pvar]
        end
        for dvar in [
            "mu_va_diff_lb",
            "mu_va_diff_ub",
            "lam_ohm_active_fr",
            "lam_ohm_active_to",
            "lam_ohm_reactive_fr",
            "lam_ohm_reactive_to",
        ]
            dres_h5[dvar][e] = brsol[dvar]
        end

        # conic duals (unrolled)
        dres_h5["nu_sm_to"][e, :] .= brsol["nu_sm_to"]
        dres_h5["nu_sm_fr"][e, :] .= brsol["nu_sm_fr"]
    end
    for p in keys(sol["buspairs"])
        bpsol = sol["buspairs"]["$p"]

        for pvar in ["wr", "wi"]
            pres_h5[pvar][p] = bpsol[pvar]
        end
    end

    dres_h5["S"] .= sol["S"]

    return res_h5
end
