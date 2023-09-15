import SparseArrays: SparseMatrixCSC, findnz

struct StandardFormData
    A::SparseMatrixCSC
    b::Vector{Float64}
    c::Vector{Float64}
    c0::Float64
    l::Vector{Float64}
    u::Vector{Float64}
    columns::Dict{VariableRef, Union{Int, Dict{String, Float64}}}
end

mutable struct StandardFormOPFModel{OPF <: PM.AbstractPowerModel}
    data::Dict{String,Any}
    model::JuMP.Model
    std::StandardFormData
    objective_kind::String
    mu::Float64
end

StandardFormDCPPowerModel = StandardFormOPFModel{PM.DCPPowerModel}

"""
    make_standard_form_data(lp::Model)

Convert a JuMP Model to standard form matrices:

    min  ∑ᵢ cᵢxᵢ + c₀
    s.t.     Ax == b
         l <= x <= u

The input model must only have linear constraints (GenericAffExpr).
If the input model has a quadratic (QuadExpr) objective, only the linear parts are used.
"""
function make_standard_form_data(lp::Model)
    jump_std = JuMP._standard_form_matrix(lp)  # NOTE: this is not an unstable public API!

    function deletecol(A, col)
        dim = length(size(A))
        if dim == 1
            A = [A[1:col-1]; A[col+1:end]]
        elseif dim == 2
            A = [A[:, 1:col-1] A[:, col+1:end]]
        else
            error("deletecol: only 1D and 2D arrays are supported")
        end
        return A
    end

    columns = Dict{VariableRef, Union{Int, Dict{String, Float64}}}()
    for (var, col) in jump_std.columns
        columns[var] = col
    end

    A = copy(jump_std.A)
    l = copy(jump_std.lower)
    u = copy(jump_std.upper)
    b = zeros(size(jump_std.A, 1))

    equal_bound_idxs = findall(l .== u)
    if length(equal_bound_idxs) > 0
        for eqb_idx in reverse(sort(equal_bound_idxs))
            x_rows = findall(A[:, eqb_idx] .!= 0)
            x_coeffs = A[x_rows, eqb_idx]
            fixed_val = l[eqb_idx]

            b[x_rows] .-=  (x_coeffs .* fixed_val)

            A = deletecol(A, eqb_idx)
            l = deletecol(l, eqb_idx)
            u = deletecol(u, eqb_idx)

            for (var, col) in columns
                if !isa(col, Int)
                    continue
                elseif col > eqb_idx
                    columns[var] = col - 1
                elseif col == eqb_idx
                    columns[var] = Dict{String, Float64}()
                    columns[var]["fixed"] = fixed_val
                end
            end
        end
    end

    n_x = size(A, 2)

    obj_sense = objective_sense(lp)
    obj_func = objective_function(lp)

    function get_c(obj_func::JuMP.AffExpr, n_x::Int)
        vars = collect(keys(obj_func.terms))
        coeffs = collect(values(obj_func.terms))

        c = zeros(n_x)
        c0 = obj_func.constant

        for (var, coeff) in zip(vars, coeffs)
            col = columns[var]
            if isa(col, Int)
                c[col] = coeff
            elseif isa(col, Dict) && haskey(col, "fixed")
                c0 += coeff * col["fixed"]
            else
                error("get_c: unexpected col")
            end
        end
        return c, c0
    end

    function get_c(obj_func::JuMP.QuadExpr, n_x::Int)
        return get_c(obj_func.aff, n_x)
    end

    c, c0 = get_c(obj_func, n_x)

    if obj_sense == JuMP.MOI.MAX_SENSE
        c = -c
        c0 = -c0
    end

    return StandardFormData(A, b, c, c0, l, u, columns)
end

"""
    make_standard_form(lp::Model, optimizer=Ipopt.Optimizer; objective_kind="linear", mu=0.1)

Given a linear JuMP model, convert it to a JuMP model in standard form:

    min  ∑ᵢ cᵢxᵢ + c₀
    s.t.     Ax == b
         l <= x <= u

The input model must only have linear constraints (`GenericAffExpr`).
If the input model has a quadratic (`QuadExpr`) objective, only the linear parts are used.

The output model can be given a barrier term using MOI.ExponentialCone.
Set `objective_kind` to "conic" to use a barrier with parameter `mu`.
"""
function make_standard_form(lp::Model, optimizer; objective_kind="linear", mu=0.1)
    std = make_standard_form_data(lp)
    N = size(std.A, 2)

    model = Model(optimizer)

    @variable(model, std.l[i] <= x[i=1:N] <= std.u[i])

    @constraint(model, constraints, std.A * x .== std.b)

    if objective_kind == "linear"
        model.ext[:objective_kind] = objective_kind
        model.ext[:mu] = 0.0

        @objective(model, Min, sum(std.c[i] * x[i] for i in 1:N) + std.c0)
    elseif objective_kind == "conic"
        model.ext[:objective_kind] = objective_kind
        model.ext[:mu] = mu

        finite_ls = isfinite.(std.l)
        finite_us = isfinite.(std.u)

        N_finite_l = sum(finite_ls)
        N_finite_u = sum(finite_us)

        finite_l_map = Dict{Int, Int}()
        finite_u_map = Dict{Int, Int}()
        finite_l_idx = 0
        finite_u_idx = 0
        for i in 1:N
            if finite_ls[i]
                finite_l_idx += 1
                finite_l_map[finite_l_idx] = i
            end
            if finite_us[i]
                finite_u_idx += 1
                finite_u_map[finite_u_idx] = i
            end
        end

        @variable(model, t_l[1:N_finite_l])
        @variable(model, t_u[1:N_finite_u])
        @constraint(model, cone_lower[i=1:N_finite_l], [t_l[i], 1, x[finite_l_map[i]] - std.l[finite_l_map[i]]] in MOI.ExponentialCone())
        @constraint(model, cone_upper[i=1:N_finite_u], [t_u[i], 1, std.u[finite_u_map[i]] - x[finite_u_map[i]]] in MOI.ExponentialCone())

        @objective(model, Min,
            sum(std.c[i] * x[i] for i in 1:N) + std.c0
            - mu * sum(t_l[i] for i in 1:N_finite_l)
            - mu * sum(t_u[i] for i in 1:N_finite_u)
        )
    else
        error("make_standard_form: unknown objective_kind")
    end

    return model, std
end

"""
    map_standard_form_solution(model::Model, columns::Dict)

Return a mapping of the original variables to their values
as given by the solution of the standard form model.
"""
function map_standard_form_solution(model::Model, columns::Dict)
    x_sol = value.(model[:x])

    solution = Dict{VariableRef, Float64}()
    for (var, col) in columns
        if isa(col, Int)
            solution[var] = x_sol[col]
        elseif isa(col, Dict) && haskey(col, "fixed")
            solution[var] = col["fixed"]
        end
    end

    return solution
end

function map_standard_form_solution(model::Model, std::StandardFormData)
    return map_standard_form_solution(model, std.columns)
end

function map_standard_form_solution(opf::StandardFormOPFModel{OPF}) where {OPF <: PM.AbstractPowerModel}
    return map_standard_form_solution(opf.model, opf.std)
end

function standard_form_data_to_dict(std::StandardFormData)
    str_columns = Dict(string(k) => v for (k,v) in std.columns)
    d = Dict{String, Any}()

    d["m"] = size(std.A, 1)
    d["n"] = size(std.A, 2)

    d["b"] = std.b
    d["c"] = std.c
    d["c0"] = std.c0

    # TODO: not able to reproduce A@x-b=0 in torch/ml4opf
    d["A_coo"] = coo = Dict{String, Any}()
    I, J, V = findnz(std.A)
    coo["I"] = I
    coo["J"] = J
    coo["V"] = V

    finite_ls = isfinite.(std.l)
    finite_us = isfinite.(std.u)

    @assert all(std.l[.!finite_ls] .== -Inf)
    @assert all(std.u[.!finite_us] .== Inf)

    d["l"] = Dict{String, Any}()
    d["l"]["finite"] = std.l[finite_ls]
    d["l"]["mask_idx"] = Int.(finite_ls)

    d["u"] = Dict{String, Any}()
    d["u"]["finite"] = std.u[finite_us]
    d["u"]["mask_idx"] = Int.(finite_us)

    d["columns"] = str_columns

    return d
end

function standard_form_data_to_dict(opf::StandardFormOPFModel{OPF}) where {OPF <: PM.AbstractPowerModel}
    return standard_form_data_to_dict(opf.std)
end

function write_standard_form_data(config::Dict, D::Dict)
    for opf_config in values(config["OPF"])
        OPF = OPFGenerator.OPF2TYPE[opf_config["type"]]
        if (OPF <: StandardFormOPFModel)
            PM.silence()
            data = make_basic_network(pglib(config["ref"]))
            opf = OPFGenerator.build_opf(OPF, data, nothing)

            Dopf = get!(D, opf_config["type"], Dict{String,Any}())
            Dopf["standard_form"] = OPFGenerator.standard_form_data_to_dict(opf)
        end
    end
end