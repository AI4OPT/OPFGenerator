"""
    build_dcopf(data, optimizer)

Build a DC-OPF model.
"""
function build_dcopf(data::Dict{String,Any}, optimizer)
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

    #
    #   I. Variables
    #

    # bus voltage angle
    JuMP.@variable(model, va[1:N])

    # generator active dispatch
    JuMP.@variable(model, ref[:gen][g]["pmin"] <= pg[g in 1:G] <= ref[:gen][g]["pmax"])

    # branch flows
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= pf[(l,i,j) in ref[:arcs_from]] <= ref[:branch][l]["rate_a"])

    # https://github.com/lanl-ansi/PowerModels.jl/blob/d9c00d228a4e875ac3426425ad6c8f8338309efa/src/form/dcp.jl#L255-L258
    # pf_expr = Dict{Any, Any}( ((l,i,j), pf[(l,i,j)]) for (l,i,j) in ref[:arcs_from] )
    # pf_expr = merge(pf_expr, Dict( ((l,j,i), -1.0*pf[(l,i,j)]) for (l,i,j) in ref[:arcs_from] ))
    pf_expr = [((l,i,j), pf[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]
    pf_expr = vcat(pf_expr, [((l,j,i), -1.0*pf[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])

    model[:pf] = JuMP.Containers.DenseAxisArray(collect(getfield.(pf_expr, 2)), collect(getfield.(pf_expr, 1)))
    pf = model[:pf]

    #
    #   II. Constraints
    #

    # Slack bus
    JuMP.@constraint(model, slack_bus, va[i0] == 0.0)

    # Nodal power balance
    JuMP.@constraint(model,
        kirchhoff[i in 1:N],
        sum(pf[a] for a in ref[:bus_arcs][i]) ==
        sum(pg[g] for g in ref[:bus_gens][i]) -
        sum(load["pd"] for load in bus_loads[i]) -
        sum(shunt["gs"] for shunt in bus_shunts[i])*1.0^2
    )

    # Branch power flow physics and limit constraints
    # We pre-allocate constraint containers for simplicity
    model[:voltage_difference_limit] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:ohm_eq] = Vector{JuMP.ConstraintRef}(undef, E)
    for (i,branch) in ref[:branch]
        data["branch"]["$i"]["br_status"] == 0 && continue  # skip branch

        f_idx = (i, branch["f_bus"], branch["t_bus"])

        p_fr = pf[f_idx]

        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PM.calc_branch_y(branch)

        # From side of the branch flow
        model[:ohm_eq][i] = JuMP.@constraint(model, p_fr == -b*(va_fr - va_to))

        # Voltage angle difference limit
        model[:voltage_difference_limit][i] = JuMP.@constraint(model, branch["angmin"] <= va_fr - va_to <= branch["angmax"])
    end

    #
    #   III. Objective
    #
    JuMP.@objective(model, Min, sum(
        gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3]
        for (i,gen) in ref[:gen]
    ))

    return model
end

"""
    _extract_dcopf_solution(model, data)

Extract DCOPF solution from optimization model.
The model must have been solved before.
"""
function _extract_dcopf_solution(model::JuMP.Model, data::Dict{String,Any})
    # Pre-process data
    ref = PM.build_ref(data)[:it][:pm][:nw][0]
    N = length(ref[:bus])
    G = length(ref[:gen])
    E = length(data["branch"])

    # Build the solution dictionary
    res = Dict{String,Any}()
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
            "va" => value(model[:va][bus]),
            # dual vars
            "lam_pb" => dual(model[:kirchhoff][bus])
        )
    end

    for b in 1:E
        if data["branch"]["$b"]["br_status"] == 0
            # branch is under outage --> we set everything to zero
            sol["branch"]["$b"] = Dict(
                "pf" => 0.0,
                # dual vars
                "mu_va_diff" => 0.0,
                "lam_ohm" => 0.0,
            )
        else
            sol["branch"]["$b"] = Dict(
                "pf" => value(model[:pf][(b,ref[:branch][b]["f_bus"],ref[:branch][b]["t_bus"])]),
                # dual vars
                "mu_va_diff" => dual(model[:voltage_difference_limit][b]),
                "lam_ohm" => dual(model[:ohm_eq][b]),
            )
        end
    end

    for g in 1:G
        sol["gen"]["$g"] = Dict(
            "pg" => value(model[:pg][g]),
            # dual vars
            "mu_pg_lb" => dual(LowerBoundRef(model[:pg][g])),
            "mu_pg_ub" => dual(UpperBoundRef(model[:pg][g]))
        )
    end

    return res
end