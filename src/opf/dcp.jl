"""
    build_dcopf(data, optimizer)

Build a DC-OPF model.
"""
function build_opf(::Type{PM.DCPPowerModel}, data::Dict{String,Any}, optimizer;
    T=Float64,    
)
    # Cleanup and pre-process data
    PM.standardize_cost_terms!(data, order=2)
    PM.calc_thermal_limits!(data)

    # Check that slack bus is unique
    ref_buses = make_ref_buses(data)
    length(ref_buses) == 1 || error("Data has $(length(ref_buses)) slack buses (expected 1).")
    i0 = first(ref_buses)

    # Grab some data
    N = length(data["bus"])
    G = length(data["gen"])
    E = length(data["branch"])

    bus_loads = make_bus_loads(data)
    bus_shunts = make_bus_shunts(data)
    bus_gens = make_bus_gens(data)

    arcs_from = make_arcs_from(data)
    bus_arcs = make_bus_arcs(data, arcs_from)

    model = JuMP.GenericModel{T}(optimizer)
    model.ext[:opf_model] = PM.DCPPowerModel  # for internal checks

    #
    #   I. Variables
    #

    # bus voltage angle
    JuMP.@variable(model, va[1:N])

    # generator active dispatch
    JuMP.@variable(model, data["gen"]["$g"]["pmin"] <= pg[g in 1:G] <= data["gen"]["$g"]["pmax"])

    # branch flows
    JuMP.@variable(model, -data["branch"]["$l"]["rate_a"] <= pf[(l,i,j) in arcs_from] <= data["branch"]["$l"]["rate_a"])

    # https://github.com/lanl-ansi/PowerModels.jl/blob/d9c00d228a4e875ac3426425ad6c8f8338309efa/src/form/dcp.jl#L255-L258
    # pf_expr = Dict{Any, Any}( ((l,i,j), pf[(l,i,j)]) for (l,i,j) in ref[:arcs_from] )
    # pf_expr = merge(pf_expr, Dict( ((l,j,i), -1.0*pf[(l,i,j)]) for (l,i,j) in ref[:arcs_from] ))
    pf_expr = [((l,i,j), pf[(l,i,j)]) for (l,i,j) in arcs_from]
    pf_expr = vcat(pf_expr, [((l,j,i), -1.0*pf[(l,i,j)]) for (l,i,j) in arcs_from])

    model[:pf] = pf = JuMP.Containers.DenseAxisArray(collect(getfield.(pf_expr, 2)), collect(getfield.(pf_expr, 1)))

    # fix disabled generators
    zero_variables([data["gen"]["$g"]["gen_status"] for g in 1:G], pg)
    # fix disabled branches
    zero_variables([data["branch"]["$(a[1])"]["br_status"] for a in arcs_from], [pf[a] for a in arcs_from])

    #
    #   II. Constraints
    #

    # Slack bus
    JuMP.@constraint(model, slack_bus, va[i0] == 0.0)

    # Nodal power balance
    JuMP.@constraint(model,
        kirchhoff[i in 1:N],
        sum(pg[g] for g in bus_gens[i] if data["gen"]["$g"]["gen_status"] == 1)
        - sum(pf[a] for a in bus_arcs[i] if data["branch"]["$(a[1])"]["br_status"] == 1)
        ==
        sum(data["load"]["$l"]["pd"] for l in bus_loads[i])
        + sum(data["shunt"]["$s"]["gs"] for s in bus_shunts[i])
    )

    # Branch power flow physics and limit constraints
    # We pre-allocate constraint containers for simplicity
    model[:voltage_difference_limit] = Vector{JuMP.ConstraintRef}(undef, E)
    model[:ohm_eq] = Vector{JuMP.ConstraintRef}(undef, E)

    for i in make_active_branches(data)
        branch = data["branch"]["$i"]

        f_idx = (i, branch["f_bus"], branch["t_bus"])

        p_fr = pf[f_idx]

        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PM.calc_branch_y(branch)

        # From side of the branch flow
        model[:ohm_eq][i] = JuMP.@constraint(model, -b*(va_fr - va_to) - p_fr == 0)

        # Voltage angle difference limit
        model[:voltage_difference_limit][i] = JuMP.@constraint(model, branch["angmin"] <= va_fr - va_to <= branch["angmax"])
    end

    #
    #   III. Objective
    #
    l, u = extrema(gen["cost"][1] for gen in values(data["gen"]) if gen["gen_status"] == 1)
    (l == u == 0.0) || @warn "Data $(data["name"]) has quadratic cost terms; those terms are being ignored"
    JuMP.@objective(model, Min, sum(
        data["gen"]["$g"]["cost"][2]*pg[g] + data["gen"]["$g"]["cost"][3]
        for g in 1:G if data["gen"]["$g"]["gen_status"] == 1
    ))

    return OPFModel{PM.DCPPowerModel}(data, model)
end

function update!(opf::OPFModel{PM.DCPPowerModel}, data::Dict{String,Any})
    opf.data = data

    bus_loads = make_bus_loads(data)
    bus_shunts = make_bus_shunts(data)
    N = length(data["bus"])

    pd = [sum(data["load"]["$l"]["pd"] for l in bus_loads[i]; init=0.0) for i in 1:N]
    gs = [sum(data["shunt"]["$s"]["gs"] for s in bus_shunts[i]; init=0.0) for i in 1:N]
    
    JuMP.set_normalized_rhs.(opf.model[:kirchhoff], pd .+ gs)

    return nothing
end

"""
    extract_result(opf::OPFModel{PM.DCPPowerModel})

Extract DCOPF solution from optimization model.
The model must have been solved before.
"""
function extract_result(opf::OPFModel{PM.DCPPowerModel})
    data  = opf.data
    model = opf.model

    # Pre-process data
    N = length(data["bus"])
    G = length(data["gen"])
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

    for bus in 1:N
        sol["bus"]["$bus"] = Dict(
            "va" => value(model[:va][bus]),
            # dual vars
            "lam_kirchhoff" => dual(model[:kirchhoff][bus])
        )
    end

    for b in 1:E
        if data["branch"]["$b"]["br_status"] == 0
            # branch is under outage --> we set everything to zero
            sol["branch"]["$b"] = Dict(
                "pf" => 0.0,
                # dual vars
                "mu_va_diff" => 0.0,
                "lam_ohm"    => 0.0,
                "mu_sm_lb"   => 0.0,
                "mu_sm_ub"   => 0.0,
            )
        else
            lij = (b, data["branch"]["$b"]["f_bus"], data["branch"]["$b"]["t_bus"])
            sol["branch"]["$b"] = Dict(
                "pf" => value(model[:pf][lij]),
                # dual vars
                "mu_va_diff" => dual(model[:voltage_difference_limit][b]),
                "lam_ohm"    => dual(model[:ohm_eq][b]),
                "mu_sm_lb"   => dual(LowerBoundRef(model[:pf][lij])),
                "mu_sm_ub"   => dual(UpperBoundRef(model[:pf][lij])),
            )
        end
    end

    for g in 1:G
        if data["gen"]["$g"]["gen_status"] == 0
            # generator is under outage --> we set everything to zero
            sol["gen"]["$g"] = Dict(
                "pg" => 0.0,
                # dual vars
                "mu_pg_lb" => 0.0,
                "mu_pg_ub" => 0.0,
            )
        else
            sol["gen"]["$g"] = Dict(
                "pg" => value(model[:pg][g]),
                # dual vars
                "mu_pg_lb" => dual(LowerBoundRef(model[:pg][g])),
                "mu_pg_ub" => dual(UpperBoundRef(model[:pg][g])),
            )
        end
    end 

    # sol["global"] = Dict(
    #     "lam_slack_bus" => dual(model[:slack_bus]),
    # )

    return res
end

function json2h5(::Type{PM.DCPPowerModel}, res)
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
        "va" => zeros(Float64, N),
        "pg" => zeros(Float64, G),
        "pf" => zeros(Float64, E),
    )
    res_h5["dual"] = dres_h5 = Dict{String,Any}(
        "lam_kirchhoff"   => zeros(Float64, N),
        "mu_pg_lb"               => zeros(Float64, G),
        "mu_pg_ub"               => zeros(Float64, G),
        "lam_ohm"                => zeros(Float64, E),
        "mu_sm_lb"               => zeros(Float64, E),
        "mu_sm_ub"               => zeros(Float64, E),
        "mu_va_diff"             => zeros(Float64, E),
        "lam_slack_bus"          => 0.0,
    )

    # extract from solution
    for i in 1:N
        bsol = sol["bus"]["$i"]

        pres_h5["va"][i] = bsol["va"]

        dres_h5["lam_kirchhoff"][i] = bsol["lam_kirchhoff"]
    end
    for g in 1:G
        gsol = sol["gen"]["$g"]

        pres_h5["pg"][g] = gsol["pg"]

        dres_h5["mu_pg_lb"][g] = gsol["mu_pg_lb"]
        dres_h5["mu_pg_ub"][g] = gsol["mu_pg_ub"]
    end
    for e in 1:E
        brsol = sol["branch"]["$e"]

        pres_h5["pf"][e] = brsol["pf"]
        
        dres_h5["lam_ohm"][e] = brsol["lam_ohm"]
        dres_h5["mu_sm_lb"][e] = brsol["mu_sm_lb"]
        dres_h5["mu_sm_ub"][e] = brsol["mu_sm_ub"]
        dres_h5["mu_va_diff"][e] = brsol["mu_va_diff"]
    end

    # dres_h5["lam_slack_bus"] = sol["global"]["lam_slack_bus"]

    return res_h5
end
