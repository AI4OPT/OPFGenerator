using JuMP
using LinearAlgebra

"""
    make_standard_form_matrix(lp::Model)

Convert a JuMP Model to standard form matrices:

    min   c ⋅ x + c0
    s.t.     Ax == b
         l <= x <= u

The input model must only have linear constraints (GenericAffExpr).
If the input model has a quadratic (QuadExpr) objective, only the linear parts are used.
"""
function make_standard_form_matrix(lp::Model)
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

    A = deepcopy(jump_std.A)
    l = deepcopy(jump_std.lower)
    u = deepcopy(jump_std.upper)
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

    return (
        A = A, b = b,
        c = c, c0 = c0,
        l = l, u = u,
        columns = columns
    )
end

function make_standard_form(opf::OPFModel{PM.DCPPowerModel}, optimizer=Ipopt.Optimizer)
    return make_standard_form(opf.model, optimizer)
end

"""
    make_standard_form(lp::Model; optimizer=Ipopt.Optimizer)

Given a linear JuMP model, convert it to a JuMP model in standard form:

    min   c ⋅ x + c0
    s.t.     Ax == b
         l <= x <= u

The input model must only have linear constraints (GenericAffExpr).
If the input model has a quadratic (QuadExpr) objective, only the linear parts are used.
"""
function make_standard_form(lp::Model, optimizer=Ipopt.Optimizer)
    std = make_standard_form_matrix(lp)

    model = Model(optimizer)

    @variable(model, std.l[i] <= x[i=1:size(std.A, 2)] <= std.u[i])

    @constraint(model, constraints, std.A * x .== std.b)

    @objective(model, Min, dot(std.c, x) + std.c0)

    return model, std.columns
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
