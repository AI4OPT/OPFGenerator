"""
    build_dcopf(data, optimizer)

Build a Standard Form model.
"""
function build_opf(::Type{StandardFormDCPPowerModel}, data::Dict{String,Any}, optimizer; kwargs...)
    opf = build_opf(PM.DCPPowerModel, data, optimizer)

    model, std = make_standard_form(opf, optimizer; kwargs...)

    model.ext[:opf_model] = PM.DCPPowerModel # NOTE: should be StandardFormDCPPowerModel?

    return StandardFormDCPPowerModel(data, model, std, model.ext[:objective_kind], model.ext[:mu])
end

"""
    _extract_dcopf_solution(model, data)

Extract standard form solution from optimization model.
The model must have been solved before.
"""
function extract_result(opf::StandardFormDCPPowerModel)
    data  = opf.data
    model = opf.model

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

    res["objective_kind"] = model.ext[:objective_kind]
    res["mu"] = model.ext[:mu]

    res["b"] = opf.std.b

    res["solution"] = sol = Dict{String,Any}()

    sol["per_unit"] = get(data, "per_unit", false)
    sol["baseMVA"]  = get(data, "baseMVA", 100.0)

    sol["x"] = value.(model[:x])
    sol["lambda"] = dual.(model[:constraints])
    
    # bound duals default to 0 -> if the lower/upper bound is ±∞, dual is fixed at 0
    # consistent with Mosek https://docs.mosek.com/latest/dotnetapi/prob-def-linear.html
    Nx = length(model[:x])
    sol["mu_lower"] = zeros(Float64, Nx)
    sol["mu_upper"] = zeros(Float64, Nx)

    for i in 1:Nx
        if has_lower_bound(model[:x][i])
            sol["mu_lower"][i] = dual(LowerBoundRef(model[:x][i]))
        end
        if has_upper_bound(model[:x][i])
            sol["mu_upper"][i] = dual(UpperBoundRef(model[:x][i]))
        end
    end

    orig_sol = OPFGenerator.map_standard_form_solution(opf)
    orig_sol = Dict(string(k) => v for (k,v) in orig_sol) # this can be better.. they have internal functions we can call for this

    sol["columns"] = Dict(string(k) => v for (k,v) in opf.std.columns)

    sol["bus"] = Dict{String,Any}()
    sol["branch"] = Dict{String,Any}()
    sol["gen"] = Dict{String,Any}()

    for bus in 1:N
        sol["bus"]["$bus"] = Dict(
            # "va" => value(model[:va][bus]),
            "va" => orig_sol["va[$bus]"],
            # TODO: dual vars
            # "lam_kirchhoff" => dual(model[:kirchhoff][bus])
        )
    end

    for b in 1:E
        if data["branch"]["$b"]["br_status"] == 0
            # branch is under outage --> we set everything to zero
            sol["branch"]["$b"] = Dict(
                "pf" => 0.0,
                # TODO: dual vars
                # "mu_va_diff" => 0.0,
                # "lam_ohm"    => 0.0,
                # "mu_sm_lb"   => 0.0,
                # "mu_sm_ub"   => 0.0,
            )
        else
            branch = ref[:branch][b]
            f_idx = (b, branch["f_bus"], branch["t_bus"])

            sol["branch"]["$b"] = Dict(
                "pf" => orig_sol["pf[$f_idx]"],

                # TODO: dual vars
                # "mu_va_diff" => dual(model[:voltage_difference_limit][b]),
                # "lam_ohm"    => dual(model[:ohm_eq][b]),
                # "mu_sm_lb"   => dual(LowerBoundRef(model[:pf][f_idx])),
                # "mu_sm_ub"   => dual(UpperBoundRef(model[:pf][f_idx])),
            )
        end
    end

    for g in 1:G
        sol["gen"]["$g"] = Dict(
            "pg" => orig_sol["pg[$g]"],

            # TODO: dual vars
            # "mu_pg_lb" => dual(LowerBoundRef(model[:pg][g])),
            # "mu_pg_ub" => dual(UpperBoundRef(model[:pg][g])),
        )
    end

    return res
end

function json2h5(::Type{StandardFormDCPPowerModel}, res)
    sol = res["solution"]
    N = length(sol["bus"])
    E = length(sol["branch"])
    G = length(sol["gen"])
    Nx = length(sol["x"])
    Nλ = length(sol["lambda"])

    res_h5 = Dict{String,Any}(
        "meta" => Dict{String,Any}(
            "termination_status" => res["termination_status"],
            "primal_status" => res["primal_status"],
            "dual_status" => res["dual_status"],
            "solve_time" => res["solve_time"],
        ),
    )

    res_h5["primal"] = pres_h5 = Dict{String,Any}(
        "x" => zeros(Float64, Nx),
         # TODO: b is not a primal variable + what if c/A change?
        "b" => zeros(Float64, Nλ),

        "va" => zeros(Float64, N),
        "pg" => zeros(Float64, G),
        "pf" => zeros(Float64, E),
    )
    res_h5["dual"] = dres_h5 = Dict{String,Any}(
        "lambda" => zeros(Float64, Nλ),
        "mu_lower" => zeros(Float64, Nx),
        "mu_upper" => zeros(Float64, Nx),

        # TODO: dual vars
        # "lam_kirchhoff"   => zeros(Float64, N),
        # "mu_pg_lb"               => zeros(Float64, G),
        # "mu_pg_ub"               => zeros(Float64, G),
        # "lam_ohm"                => zeros(Float64, E),
        # "mu_sm_lb"               => zeros(Float64, E),
        # "mu_sm_ub"               => zeros(Float64, E),
        # "mu_va_diff"             => zeros(Float64, E),
    )

    # extract from solution
    for i in 1:N
        bsol = sol["bus"]["$i"]

        pres_h5["va"][i] = bsol["va"]
        
        # TODO: dual vars
        # dres_h5["lam_kirchhoff"][i] = bsol["lam_kirchhoff"]
    end
    for g in 1:G
        gsol = sol["gen"]["$g"]

        pres_h5["pg"][g] = gsol["pg"]
        
        # TODO: dual vars
        # dres_h5["mu_pg_lb"][g] = gsol["mu_pg_lb"]
        # dres_h5["mu_pg_ub"][g] = gsol["mu_pg_ub"]
    end
    for e in 1:E
        brsol = sol["branch"]["$e"]

        pres_h5["pf"][e] = brsol["pf"]
        
        # TODO: dual vars
        # dres_h5["lam_ohm"][e] = brsol["lam_ohm"]
        # dres_h5["mu_sm_lb"][e] = brsol["mu_sm_lb"]
        # dres_h5["mu_sm_ub"][e] = brsol["mu_sm_ub"]
        # dres_h5["mu_va_diff"][e] = brsol["mu_va_diff"]
    end
    for x in 1:Nx
        pres_h5["x"][x] = sol["x"][x]
        
        dres_h5["mu_lower"][x] = sol["mu_lower"][x]
        dres_h5["mu_upper"][x] = sol["mu_upper"][x]
    end
    for λ in 1:Nλ
        pres_h5["b"][λ] = res["b"][λ]

        dres_h5["lambda"][λ] = sol["lambda"][λ]
    end

    return res_h5
end
