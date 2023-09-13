import SparseArrays: SparseMatrixCSC

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
end

StandardFormDCPPowerModel = StandardFormOPFModel{PM.DCPPowerModel}

"""
    make_standard_form_data(lp::Model)

Convert a JuMP Model to standard form matrices:

    min   c ⋅ x + c0
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

    return StandardFormData(A, b, c, c0, l, u, columns)
end

"""
    make_standard_form(lp::Model, optimizer=Ipopt.Optimizer; objective_type="linear", mu=0.1)

Given a linear JuMP model, convert it to a JuMP model in standard form:

    min   c ⋅ x + c0
    s.t.     Ax == b
         l <= x <= u

The input model must only have linear constraints (`GenericAffExpr`).
If the input model has a quadratic (`QuadExpr`) objective, only the linear parts are used.

The output model can be given a barrier term using an ExponentialCone.
Set `objective_type` to "cone" to use these barrier terms with parameter `mu`.
Make sure the solver supports ExponentialCone.
"""
function make_standard_form(lp::Model, optimizer; objective_type="linear", mu=0.1)
    std = make_standard_form_data(lp)
    N = size(std.A, 2)

    model = Model(optimizer)

    @variable(model, std.l[i] <= x[i=1:N] <= std.u[i])

    @constraint(model, constraints, std.A * x .== std.b)

    if objective_type == "linear"
        @objective(model, Min,
            sum(
                std.c[i]*x[i]
                for i in 1:N
            ) + std.c0
        )
    elseif objective_type == "cone"
        finite_u = isfinite.(std.u)
        model[:t_u] = Vector{JuMP.VariableRef}(undef, N)
        model[:cone_upper] = Vector{JuMP.ConstraintRef}(undef, N)

        finite_l = isfinite.(std.l)
        model[:t_l] = Vector{JuMP.VariableRef}(undef, N)
        model[:cone_lower] = Vector{JuMP.ConstraintRef}(undef, N)

        for i in 1:N
            if finite_u[i]
                model[:t_u][i] = @variable(model)
                model[:cone_upper][i] = @constraint(model, [model[:t_u][i], 1, std.u[i] - x[i]] in MOI.ExponentialCone())
            end
            if finite_l[i]
                model[:t_l][i] = @variable(model)
                model[:cone_lower][i] = @constraint(model, [model[:t_l][i], 1, x[i] - std.l[i]] in MOI.ExponentialCone())
            end
        end

        JuMP.@objective(model, Min,
            sum(
                std.c[i]*x[i]
                for i in 1:N
            ) + std.c0
            - mu * sum(model[:t_l][i] for i in 1:N if finite_l[i])
            - mu * sum(model[:t_u][i] for i in 1:N if finite_u[i])
        )
    else
        error("make_standard_form: unknown objective_type")
    end

    return model, std
end

function make_standard_form(opf::OPFModel{OPF}, optimizer; kwargs...) where {OPF <: PM.AbstractPowerModel}
    return make_standard_form(opf.model, optimizer; kwargs...)
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

    d["A"] = Array(std.A)
    d["b"] = std.b
    d["c"] = std.c
    d["c0"] = std.c0
    d["l"] = std.l
    d["u"] = std.u

    cols = d["columns"] = Dict{String, Any}()
    cols["keys"] = Array(collect(keys(str_columns)))
    cols["values"] = Array(collect(values(str_columns)))

    return d
end

function standard_form_data_to_dict(opf::StandardFormOPFModel{OPF}) where {OPF <: PM.AbstractPowerModel}
    return standard_form_data_to_dict(opf.std)
end

function standard_form_data_to_h5(opf::StandardFormOPFModel{OPF}, filename::AbstractString) where {OPF <: PM.AbstractPowerModel}
    h5open(filename, "w") do file
        d = standard_form_data_to_dict(opf)
        gr = create_group(file, "standard_form")
        for (k, v) in d
            if k == "columns"
                gr_ = create_group(gr, k)
                gr_["keys"] = v["keys"]
                gr_["values"] = v["values"]
            else
                gr[k] = v
            end
        end
    end
end

function write_standard_form_data(config::Dict)
    for opf_str in keys(config["OPF"])
        OPF = OPFGenerator.OPF2TYPE[opf_str]
        if (OPF <: StandardFormOPFModel)
            data = make_basic_network(pglib(config["ref"]))
            opf = OPFGenerator.build_opf(OPF, data, nothing)
            OPFGenerator.standard_form_data_to_h5(opf, joinpath(config["export_dir"], "$opf_str.h5"))
        end
    end
end