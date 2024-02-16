using SparseArrays

mutable struct DCPPTDFPowerModel <: PM.AbstractDCPModel end

"""
    build_dcopf(data, optimizer)

Build a DC-OPF model using PTDF (no reserves).
"""
function build_opf(::Type{DCPPTDFPowerModel}, data::Dict{String,Any}, optimizer)
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

    model = JuMP.Model(optimizer)
    model.ext[:opf_model] = DCPPTDFPowerModel  # for internal checks
    model.ext[:formulation] = "PTDF"
    model.ext[:PTDF] = PM.calc_basic_ptdf_matrix(data)

    #
    #   I. Variables
    #

    # generator active dispatch
    JuMP.@variable(model, ref[:gen][g]["pmin"] <= pg[g in 1:G] <= ref[:gen][g]["pmax"])
    JuMP.@variable(model, -ref[:branch][e]["rate_a"] <= pf[e in 1:E] <= ref[:branch][e]["rate_a"])

    JuMP.@constraint(model,
        power_balance,
        sum(pg[g] for g in 1:G) == sum(ref[:load][l]["pd"] for l in 1:length(ref[:load]))
    )

    # to be added during solve!
    model[:ptdf_flow] = Vector{JuMP.ConstraintRef}(undef, E)
    model.ext[:tracked_branches] = zeros(Bool, E)

    #
    #   III. Objective
    #

    l, u = extrema(gen["cost"][1] for (i, gen) in ref[:gen])
    (l == u == 0.0) || @warn "Data $(data["name"]) has quadratic cost terms; those terms are being ignored"

    l, u = extrema(gen["cost"][3] for (i, gen) in ref[:gen])
    (l == u == 0.0) || @warn "Data $(data["name"]) has constant cost terms; those terms are being ignored"

    JuMP.@objective(model, Min, sum(
        gen["cost"][2] * pg[i]
        for (i, gen) in ref[:gen]
    ))

    return OPFModel{DCPPTDFPowerModel}(data, model)
end


function solve!(opf::OPFModel{DCPPTDFPowerModel})
    if !haskey(opf.model.ext, :formulation) || opf.model.ext[:formulation] != "PTDF"
        error("Expected PTDF formulation")
    end

    data = opf.data
    model = opf.model

    N = length(data["bus"])
    G = length(data["gen"])
    E = length(data["branch"])
    L = length(data["load"])
    Al = sparse(
        [data["load"]["$l"]["load_bus"] for l in 1:L],
        [l for l in 1:L],
        ones(Float64, L),
        N,
        L,
    )

    Ag = sparse(
        [data["gen"]["$g"]["gen_bus"] for g in 1:G],
        [g for g in 1:G],
        ones(Float64, G),
        N,
        G,
    )

    pd = [data["load"]["$l"]["pd"] for l in 1:L]
    rate_a = [data["branch"]["$e"]["rate_a"] for e in 1:E]

    p_expr = Ag * model[:pg] - Al * pd

    solved = false
    niter = 0
    while !solved
        optimize!(opf.model, _differentiation_backend = MathOptSymbolicAD.DefaultBackend())
        st = termination_status(model)
        st == MOI.OPTIMAL || error("Solver failed to converge: $st")
        pg_ = value.(model[:pg])
        p_ = Ag * pg_ - Al * pd
        pf_ = model.ext[:PTDF] * p_

        n_violated = 0
        for e in 1:E
            if (model.ext[:tracked_branches][e]) || (-rate_a[e] <= pf_[e] <= rate_a[e])
                continue
            end

            model[:ptdf_flow][e] = JuMP.@constraint(
                model,
                model[:pf][e] == dot(model.ext[:PTDF][e, :], p_expr)
            )
            model.ext[:tracked_branches][e] = true
            n_violated += 1
        end
        solved = n_violated == 0
        niter += 1
    end

    @info("Solved in $niter PTDF iteration" * (niter == 1 ? "" : "s"))

    opf.model.ext[:ptdf_iter] = niter
end


"""
    _extract_dcopf_solution(model, data)

Extract DCOPF solution from optimization model.
The model must have been solved before.
"""
function extract_result(opf::OPFModel{DCPPTDFPowerModel})
    data = opf.data
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
    sol["baseMVA"] = get(data, "baseMVA", 100.0)

    ### populate branches, gens, buses ###

    sol["branch"] = Dict{String,Any}()
    sol["gen"] = Dict{String,Any}()

    sol["power_balance"] = dual(model[:power_balance])

    for b in 1:E
        if data["branch"]["$b"]["br_status"] == 0
            # branch is under outage --> we set everything to zero
            sol["branch"]["$b"] = Dict(
                "pf" => 0.0,
                # dual vars
                "mu_pf_lb" => 0.0,
                "mu_pf_ub" => 0.0,
                "lam_ptdf_flow" => 0.0,
            )
        else
            sol["branch"]["$b"] = Dict(
                "pf" => value(model[:pf][b]),
                # dual vars
                "mu_pf_lb" => dual(LowerBoundRef(model[:pf][b])),
                "mu_pf_ub" => dual(UpperBoundRef(model[:pf][b])),
                # if not tracked, dual is 0
                "lam_ptdf_flow" => isdefined(model[:ptdf_flow], b) ? dual(model[:ptdf_flow][b]) : 0.0,
            )
        end
    end

    for g in 1:G
        sol["gen"]["$g"] = Dict(
            "pg" => value(model[:pg][g]),
            # dual vars
            "mu_pg_lb" => dual(LowerBoundRef(model[:pg][g])),
            "mu_pg_ub" => dual(UpperBoundRef(model[:pg][g])),
        )
    end

    return res
end

function json2h5(::Type{DCPPTDFPowerModel}, res)
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
        "pg" => zeros(Float64, G),
        "pf" => zeros(Float64, E),
    )
    res_h5["dual"] = dres_h5 = Dict{String,Any}(
        "mu_pg_lb" => zeros(Float64, G),
        "mu_pg_ub" => zeros(Float64, G),
        "mu_pf_lb" => zeros(Float64, E),
        "mu_pf_ub" => zeros(Float64, E),
        "lam_ptdf_flow" => zeros(Float64, E),
    )

    # extract from solution
    for g in 1:G
        gsol = sol["gen"]["$g"]

        pres_h5["pg"][g] = gsol["pg"]

        dres_h5["mu_pg_lb"][g] = gsol["mu_pg_lb"]
        dres_h5["mu_pg_ub"][g] = gsol["mu_pg_ub"]
    end
    for e in 1:E
        brsol = sol["branch"]["$e"]

        pres_h5["pf"][e] = brsol["pf"]

        dres_h5["mu_pf_lb"][e] = brsol["mu_pf_lb"]
        dres_h5["mu_pf_ub"][e] = brsol["mu_pf_ub"]
        dres_h5["lam_ptdf_flow"][e] = brsol["lam_ptdf_flow"]

    end

    dres_h5["lam_power_balance"] = sol["power_balance"]

    return res_h5
end
