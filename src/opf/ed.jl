using SparseArrays

mutable struct EconomicDispatch <: PM.AbstractPowerModel end

const POWER_BALANCE_PENALTY = 350000.0
const RESERVE_SHORTAGE_PENALTY = 110000.0
const THERMAL_PENALTY = 150000.0
const MAX_PTDF_ITERATIONS = 128
const MAX_PTDF_PER_ITERATION = 8
const ITERATIVE_PTDF_TOL = 1e-6

function build_opf(::Type{EconomicDispatch}, data::Dict{String,Any}, optimizer;
    T=Float64,
    soft_power_balance::Bool=false,
    soft_reserve_requirement::Bool=false,
    soft_thermal_limit::Bool=false,
    power_balance_penalty=POWER_BALANCE_PENALTY,
    reserve_shortage_penalty=RESERVE_SHORTAGE_PENALTY,
    thermal_penalty=THERMAL_PENALTY,
    iterative_ptdf::Bool=true,
    iterative_ptdf_tol=ITERATIVE_PTDF_TOL,
    max_ptdf_iterations=MAX_PTDF_ITERATIONS,
    max_ptdf_per_iteration=MAX_PTDF_PER_ITERATION,
)
    # Cleanup and pre-process data
    PM.standardize_cost_terms!(data, order=2)
    PM.calc_thermal_limits!(data)
    ref = PM.build_ref(data)[:it][:pm][:nw][0]

    # Check that slack bus is unique
    length(ref[:ref_buses]) == 1 || error("Data has $(length(ref[:ref_buses])) slack buses (expected 1).")
    i0 = first(keys(ref[:ref_buses]))

    power_balance_penalty >= 0.0 || error("EconomicDispatch option power_balance_penalty must be non-negative")
    reserve_shortage_penalty >= 0.0 || error("EconomicDispatch option reserve_shortage_penalty must be non-negative")
    thermal_penalty >= 0.0 || error("EconomicDispatch option transmission_penalty must be non-negative")
    max_ptdf_iterations > 0 || error("EconomicDispatch option max_ptdf_iterations must be a positive integer")
    max_ptdf_per_iteration > 0 || error("EconomicDispatch option max_ptdf_per_iteration must be a positive integer")

    # Grab some data
    N = length(ref[:bus])
    G = length(ref[:gen])
    E = length(ref[:branch])
    L = length(ref[:load])

    pmin = [ref[:gen][g]["pmin"] for g in 1:G]
    pmax = [ref[:gen][g]["pmax"] for g in 1:G]
    PD = sum(ref[:load][l]["pd"] for l in 1:L)
    pfmax = [data["branch"]["$e"]["rate_a"] for e in 1:E]

    minimum_reserve = get(data, "minimum_reserve", 0.0)
    # if minimum reserve is zero, rmin = rmax = 0 to fix the reserve variables to zero.
    # otherwise, rmin and rmax must be present in the data (not standard in PGLib)
    rmin = (minimum_reserve > 0.0) ? [ref[:gen][g]["rmin"] for g in 1:G] : zeros(T, G)
    rmax = (minimum_reserve > 0.0) ? [ref[:gen][g]["rmax"] for g in 1:G] : zeros(T, G)

    model = JuMP.GenericModel{T}(optimizer)
    model.ext[:opf_model] = EconomicDispatch  # for internal checks
    model.ext[:solve_metadata] = Dict{Symbol,Any}(
        :iterative_ptdf => iterative_ptdf,
        :iterative_ptdf_tol => iterative_ptdf_tol,
        :max_ptdf_iterations => max_ptdf_iterations,
        :max_ptdf_per_iteration => max_ptdf_per_iteration,
    )

    # Variables

    JuMP.@variable(model, pmin[g] <= pg[g in 1:G] <= pmax[g])
    JuMP.@variable(model, rmin[g] <= r[g in 1:G] <= rmax[g])
    JuMP.@variable(model, pf[e in 1:E])

    JuMP.@variable(model, δpb_surplus >= 0)
    JuMP.@variable(model, δpb_shortage >= 0)
    JuMP.@variable(model, δr_shortage >= 0)
    JuMP.@variable(model, δf[1:E] >= 0)

    JuMP.@constraint(model,
        pf_lower_bound[e in 1:E],
        pf[e] + δf[e] >= -pfmax[e] 
    )
    JuMP.@constraint(model,
        pf_upper_bound[e in 1:E],
        pf[e] - δf[e] <= pfmax[e] 
    )

    if !soft_power_balance
        JuMP.set_upper_bound(δpb_surplus, 0)
        JuMP.set_upper_bound(δpb_shortage, 0)
    end

    if !soft_reserve_requirement && minimum_reserve > 0.0
        JuMP.set_upper_bound(δr_shortage, 0)
    end

    if !soft_thermal_limit
        JuMP.set_upper_bound.(δf, 0)
    end

    # Constraints

    JuMP.@constraint(model,
        total_generation[g in 1:G],
        pg[g] + r[g] <= pmax[g]
    )

    JuMP.@constraint(model,
        power_balance,
        sum(pg) + δpb_surplus - δpb_shortage == PD
    )

    JuMP.@constraint(model,
        reserve_requirement,
        sum(r) + δr_shortage >= minimum_reserve
    )

    model.ext[:PTDF] = PM.calc_basic_ptdf_matrix(data)
    if iterative_ptdf
        model[:ptdf_flow] = Vector{JuMP.ConstraintRef}(undef, E)
        model.ext[:tracked_branches] = zeros(Bool, E)
        model.ext[:ptdf_iterations] = 0
    else
        Ag, Al, pd = _ptdf_terms_from_data(data; T=T)

        JuMP.@constraint(model, ptdf_flow, model.ext[:PTDF] * (Ag * model[:pg]) - model[:pf] .== model.ext[:PTDF] * (Al * pd))
        model.ext[:tracked_branches] = ones(Bool, E)
        model.ext[:ptdf_iterations] = -1

        (E > 1024) && @warn "EconomicDispatch: Iterative PTDF is disabled; adding all $E constraints at once. Consider using the iterative PTDF option for faster convergence on large cases."
    end

    # Objective

    costs = [ref[:gen][g]["cost"] for g in 1:G]
    l, u = extrema(costs[1] for (i, gen) in ref[:gen])
    (l == u == 0.0) || @warn "Data $(data["name"]) has quadratic cost terms; those terms are being ignored"

    JuMP.@objective(model, Min,
        sum(costs[g][2] * pg[g] + costs[g][3] for g in 1:G) +
        power_balance_penalty * δpb_surplus +
        power_balance_penalty * δpb_shortage +
        reserve_shortage_penalty * δr_shortage +
        thermal_penalty * sum(δf)
    )

    return OPFModel{EconomicDispatch}(data, model)
end

function update!(opf::OPFModel{EconomicDispatch}, data::Dict{String,Any})
    PM.standardize_cost_terms!(data, order=2)
    PM.calc_thermal_limits!(data)
    ref = PM.build_ref(data)[:it][:pm][:nw][0]

    opf.data = data
    T = typeof(opf.model).parameters[1]

    L = length(ref[:load])
    PD = sum(ref[:load][l]["pd"] for l in 1:L)

    JuMP.set_normalized_rhs(opf.model[:power_balance], PD)

    MRR = get(data, "minimum_reserve", 0.0)
    JuMP.set_normalized_rhs(opf.model[:reserve_requirement], MRR)

    G = length(ref[:gen])
    rmin = (MRR > 0.0) ? [ref[:gen][g]["rmin"] for g in 1:G] : zeros(T, G)
    rmax = (MRR > 0.0) ? [ref[:gen][g]["rmax"] for g in 1:G] : zeros(T, G)

    JuMP.set_lower_bound.(opf.model[:r], rmin)
    JuMP.set_upper_bound.(opf.model[:r], rmax)

    if opf.model.ext[:solve_metadata][:iterative_ptdf]
        opf.model.ext[:tracked_branches] .= false
        opf.model.ext[:ptdf_iterations] = 0
        
        E = length(ref[:branch])
        opf.model[:ptdf_flow] = Vector{JuMP.ConstraintRef}(undef, E)

        opf.model.ext[:termination_info] = Dict{Symbol,Any}(
            :termination_status => nothing,
            :primal_status => nothing,
            :dual_status => nothing,
            :solve_time => nothing,
        )
    else
        Ag, Al, pd = _ptdf_terms_from_data(data; T=T)

        JuMP.set_normalized_rhs.(opf.model[:ptdf_flow], PTDF * Al * pd)
    end
    
    return nothing
end

function _ptdf_terms_from_data(data::Dict{String,Any}; T=Float64)
    N = length(data["bus"])
    G = length(data["gen"])
    L = length(data["load"])
    E = length(data["branch"])

    Al = sparse(
        [data["load"]["$l"]["load_bus"] for l in 1:L],
        [l for l in 1:L],
        ones(T, L),
        N,
        L,
    )

    Ag = sparse(
        [data["gen"]["$g"]["gen_bus"] for g in 1:G],
        [g for g in 1:G],
        ones(T, G),
        N,
        G,
    )

    pd = [data["load"]["$l"]["pd"] for l in 1:L]

    return Ag, Al, pd
end

function solve!(opf::OPFModel{EconomicDispatch})
    data = opf.data
    model = opf.model
    T = typeof(model).parameters[1]
    tol = model.ext[:solve_metadata][:iterative_ptdf_tol]

    E = length(data["branch"])
    rate_a = [data["branch"]["$e"]["rate_a"] for e in 1:E]
    
    Ag, Al, pd = _ptdf_terms_from_data(data, T=T)

    solved = false
    niter = 0
    solve_time = 0.0
    model.ext[:termination_info] = Dict{Symbol,Any}(
        :termination_status => nothing,
        :primal_status => nothing,
        :dual_status => nothing,
        :solve_time => nothing,
    )
    st = nothing
    solve_time = @elapsed while !solved && niter < model.ext[:solve_metadata][:max_ptdf_iterations]
        optimize!(opf.model, _differentiation_backend = MathOptSymbolicAD.DefaultBackend())
        
        st = termination_status(model)
        st ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) || break

        pg_ = value.(model[:pg])
        p_ = Ag * pg_ - Al * pd
        pf_ = model.ext[:PTDF] * p_

        n_violated = 0
        for e in 1:E
            if (model.ext[:tracked_branches][e]) || (-rate_a[e] - tol <= pf_[e] <= rate_a[e] + tol)
                continue
            end

            n_violated += 1
            if n_violated <= model.ext[:solve_metadata][:max_ptdf_per_iteration]
                model[:ptdf_flow][e] = JuMP.@constraint(
                    model,
                    dot(model.ext[:PTDF][e, :], (Ag * model[:pg])) - model[:pf][e]
                    == dot(model.ext[:PTDF][e, :], (Al * pd))
                )
                model.ext[:tracked_branches][e] = true
            end
        end
        solved = n_violated == 0
        niter += 1
    end

    if niter == model.ext[:solve_metadata][:max_ptdf_iterations]
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.INFEASIBLE_POINT,
            :dual_status => MOI.FEASIBLE_POINT,
        )
        st = MOI.ITERATION_LIMIT
    elseif st == MOI.TIME_LIMIT
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.INFEASIBLE_POINT,
            :dual_status => MOI.FEASIBLE_POINT,
        )
    elseif st ∈ (MOI.INFEASIBLE, MOI.LOCALLY_INFEASIBLE)
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.UNKNOWN_RESULT_STATUS,
            :dual_status => MOI.INFEASIBILITY_CERTIFICATE,
        )
    elseif st == MOI.DUAL_INFEASIBLE
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.INFEASIBILITY_CERTIFICATE,
            :dual_status => MOI.UNKNOWN_RESULT_STATUS,
        )
    elseif st ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.FEASIBLE_POINT,
            :dual_status => MOI.FEASIBLE_POINT,
        )
    else
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.UNKNOWN_RESULT_STATUS,
            :dual_status => MOI.UNKNOWN_RESULT_STATUS,
        )
    end

    model.ext[:termination_info][:solve_time] = solve_time
    model.ext[:termination_info][:ptdf_iterations] = niter
    model.ext[:termination_info][:termination_status] = st

    return
end

function extract_result(opf::OPFModel{EconomicDispatch})
    data = opf.data
    model = opf.model

    ref = PM.build_ref(data)[:it][:pm][:nw][0]
    N = length(ref[:bus])
    G = length(ref[:gen])
    E = length(data["branch"])

    res = Dict{String,Any}()
    res["opf_model"] = string(model.ext[:opf_model])
    res["objective"] = JuMP.objective_value(model)
    res["objective_lb"] = -Inf
    res["optimizer"] = JuMP.solver_name(model)

    res["solve_metadata"] = model.ext[:solve_metadata]
    if model.ext[:solve_metadata][:iterative_ptdf]
        tinfo = model.ext[:termination_info]
        res["termination_status"] = tinfo[:termination_status]
        res["primal_status"] = tinfo[:primal_status]
        res["dual_status"] = tinfo[:dual_status]
        res["solve_time"] = tinfo[:solve_time]
        res["ptdf_iterations"] = tinfo[:ptdf_iterations]
    else
        res["termination_status"] = JuMP.termination_status(model)
        res["primal_status"] = JuMP.primal_status(model)
        res["dual_status"] = JuMP.dual_status(model)
        res["solve_time"] = JuMP.solve_time(model)
    end
    res["minimum_reserve"] = get(data, "minimum_reserve", 0.0)

    res["solution"] = sol = Dict{String,Any}()

    sol["per_unit"] = get(data, "per_unit", false)
    sol["baseMVA"]  = get(data, "baseMVA", 100.0)

    sol["gen"] = Dict{String,Any}()
    for g in 1:G
        sol["gen"]["$g"] = Dict(
            "pg" => value(model[:pg][g]),
            "r" => value(model[:r][g]),
            "rmin" => get(data["gen"]["$g"], "rmin", 0.0),
            "rmax" => get(data["gen"]["$g"], "rmax", 0.0),

            "mu_pg" => dual(LowerBoundRef(model[:pg][g])) - dual(UpperBoundRef(model[:pg][g])),
            "mu_r" => dual(LowerBoundRef(model[:r][g])) - dual(UpperBoundRef(model[:r][g])),

            "mu_total_generation" => dual(model[:total_generation][g]),
        )
    end

    sol["branch"] = Dict{String,Any}()
    for e in 1:E
        if data["branch"]["$e"]["br_status"] == 0
            sol["branch"]["$e"] = Dict(
                "pf" => 0,
                "df" => 0,
                "mu_pf" => 0,
                "mu_df" => 0,
                "lam_ptdf" => 0,
            )
        else
            sol["branch"]["$e"] = Dict(
                "pf" => value(model[:pf][e]),
                "df" => value(model[:δf][e]),
                "mu_pf" => dual(model[:pf_lower_bound][e]) - dual(model[:pf_upper_bound][e]),
                "mu_df" => dual(LowerBoundRef(model[:δf][e])) - (has_upper_bound(model[:δf][e]) ? dual(UpperBoundRef(model[:δf][e])) : 0.0),
                "lam_ptdf" => isdefined(model[:ptdf_flow], e) ? dual(model[:ptdf_flow][e]) : 0.0,
            )
        end
    end

    sol["singleton"] = Dict{String,Any}(
        "dpb_surplus" => value(model[:δpb_surplus]),
        "dpb_shortage" => value(model[:δpb_shortage]),
        "dr_shortage" => value(model[:δr_shortage]),

        "mu_dpb_surplus" => dual(LowerBoundRef(model[:δpb_surplus])) - (has_upper_bound(model[:δpb_surplus]) ? dual(UpperBoundRef(model[:δpb_surplus])) : 0.0),
        "mu_dpb_shortage" => dual(LowerBoundRef(model[:δpb_shortage])) - (has_upper_bound(model[:δpb_shortage]) ? dual(UpperBoundRef(model[:δpb_shortage])) : 0.0),
        "mu_dr_shortage" => dual(LowerBoundRef(model[:δr_shortage])) - (has_upper_bound(model[:δr_shortage]) ? dual(UpperBoundRef(model[:δr_shortage])) : 0.0),

        "mu_power_balance" => dual(model[:power_balance]),
        "mu_reserve_requirement" => dual(model[:reserve_requirement]),
    )

    return res
end

function json2h5(::Type{EconomicDispatch}, res)
    sol = res["solution"]

    E = length(sol["branch"])
    G = length(sol["gen"])

    res_h5 = Dict{String,Any}(
        "meta" => Dict{String,Any}(
            "termination_status" => string(res["termination_status"]),
            "primal_status" => string(res["primal_status"]),
            "dual_status" => string(res["dual_status"]),
            "solve_time" => res["solve_time"],
        ),
    )

    res_h5["primal"] = pres_h5 = Dict{String,Any}(
        "pg" => zeros(Float64, G),
        "r" => zeros(Float64, G),
        "pf" => zeros(Float64, E),
        "rmin" => zeros(Float64, G),
        "rmax" => zeros(Float64, G),

        "df" => zeros(Float64, E),
        "dpb_surplus" => 0,
        "dpb_shortage" => 0,
        "dr_shortage" => 0,
    )

    res_h5["dual"] = dres_h5 = Dict{String,Any}(
        "mu_pg" => zeros(Float64, G),
        "mu_r" => zeros(Float64, G),
        "mu_total_generation" => zeros(Float64, G),
        "mu_pf" => zeros(Float64, E),
        "mu_df" => zeros(Float64, E),
        "lam_ptdf" => zeros(Float64, E),
        "mu_dpb_surplus" => 0,
        "mu_dpb_shortage" => 0,
        "mu_dr_shortage" => 0,
        "mu_power_balance" => 0,
        "mu_reserve_requirement" => 0,
    )

    for g in 1:G
        gsol = sol["gen"]["$g"]

        pres_h5["pg"][g] = gsol["pg"]
        pres_h5["r"][g] = gsol["r"]
        pres_h5["rmin"][g] = gsol["rmin"]
        pres_h5["rmax"][g] = gsol["rmax"]

        dres_h5["mu_pg"][g] = gsol["mu_pg"]
        dres_h5["mu_r"][g] = gsol["mu_r"]
        dres_h5["mu_total_generation"][g] = gsol["mu_total_generation"]
    end

    for e in 1:E
        brsol = sol["branch"]["$e"]

        pres_h5["pf"][e] = brsol["pf"]
        pres_h5["df"][e] = brsol["df"]

        dres_h5["mu_pf"][e] = brsol["mu_pf"]
        dres_h5["mu_df"][e] = brsol["mu_df"]
        dres_h5["lam_ptdf"][e] = brsol["lam_ptdf"]
    end

    pres_h5["dpb_surplus"] = sol["singleton"]["dpb_surplus"]
    pres_h5["dpb_shortage"] = sol["singleton"]["dpb_shortage"]
    pres_h5["dr_shortage"] = sol["singleton"]["dr_shortage"]
    pres_h5["minimum_reserve"] = res["minimum_reserve"]

    dres_h5["mu_dpb_surplus"] = sol["singleton"]["mu_dpb_surplus"]
    dres_h5["mu_dpb_shortage"] = sol["singleton"]["mu_dpb_shortage"]
    dres_h5["mu_dr_shortage"] = sol["singleton"]["mu_dr_shortage"]
    dres_h5["mu_power_balance"] = sol["singleton"]["mu_power_balance"]
    dres_h5["mu_reserve_requirement"] = sol["singleton"]["mu_reserve_requirement"]

    return res_h5

end

function export_ptdf(opf::OPFModel{EconomicDispatch}, filepath)
    data = opf.data
    model = opf.model

    PTDF = model.ext[:PTDF]

    h5open(filepath, "w") do file
        write(file, "PTDF", PTDF)
    end
    return nothing
end

function export_ptdf(data::Dict{String,Any}, filepath)
    PTDF = PM.calc_basic_ptdf_matrix(data)

    h5open(filepath, "w") do file
        write(file, "PTDF", PTDF)
    end
    return nothing
end