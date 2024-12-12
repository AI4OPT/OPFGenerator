function compute_voltage_phasor_bounds(data::OPFData)
    E = data.E
    vmin, vmax = data.vmin, data.vmax
    dvamin, dvamax = data.dvamin, data.dvamax
    bus_fr, bus_to = data.bus_fr, data.bus_to

    wr_min = zeros(Float64, E)
    wr_max = zeros(Float64, E)
    wi_min = zeros(Float64, E)
    wi_max = zeros(Float64, E)
    for e in 1:E
        _wrmin, _wrmax, _wimin, _wimax = compute_voltage_phasor_bounds(
            vmin[bus_fr[e]],
            vmax[bus_fr[e]],
            vmin[bus_to[e]],
            vmax[bus_to[e]],
            dvamin[e],
            dvamax[e],
        )
        wr_min[e] = _wrmin
        wr_max[e] = _wrmax
        wi_min[e] = _wimin
        wi_max[e] = _wimax
    end

    return (wr_min, wr_max, wi_min, wi_max)
end

"""
    compute_voltage_phasor_bounds(vfmin, vfmax, vtmin, vtmax, dvamin, dvamax)

Compute lower/upper bounds on wr/wi variables.
"""
function compute_voltage_phasor_bounds(vfmin, vfmax, vtmin, vtmax, dvamin, dvamax)
    abs(dvamin) <= Base.Math.pi/2 || throw(DomainError("dvamin must be in [-π/2, π/2]"))
    abs(dvamax) <= Base.Math.pi/2 || throw(DomainError("dvamax must be in [-π/2, π/2]"))

    sinmin, cosmin = sincos(dvamin)
    sinmax, cosmax = sincos(dvamax)

    wmin = vfmin * vtmin
    wmax = vfmax * vtmax

    if dvamin >= 0
        wr_min = wmin * cosmax
        wr_max = wmax * cosmin
        wi_min = wmin * sinmin
        wi_max = wmax * sinmax
    elseif dvamax <= 0
        wr_min = wmin * cosmin
        wr_max = wmax * cosmax
        wi_min = wmax * sinmin
        wi_max = wmin * sinmax
    else
        # dvamin < 0 < dvamax
        wr_min = wmin * min(cosmin, cosmax)
        wr_max = wmax * 1.0  # max wr attained at θ = 0
        wi_min = wmax * sinmin
        wi_max = wmax * sinmax
    end

    return wr_min, wr_max, wi_min, wi_max
end

function extract_metadata(opf::OPFModel{<:AbstractFormulation})
    model = opf.model

    tst = termination_status(model)
    pst = primal_status(model)
    dst = dual_status(model)

    # Query primal/dual objective value only if 
    # * a primal/dual solution exists
    # * it's a point (i.e. not a ray, not a reduction certificate)
    OK_STATUSES = [MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT, MOI.INFEASIBLE_POINT]

    z_primal = if has_values(model) && pst ∈ OK_STATUSES
        objective_value(model)
    else
        NaN
    end
    z_dual = if has_duals(model) && dst ∈ OK_STATUSES
        dual_objective_value(model) 
    else
        NaN
    end

    return Dict{String,Any}(
        "formulation" => string(model.ext[:opf_model]),
        "termination_status" => string(tst),
        "primal_status" => string(pst),
        "dual_status" => string(dst),
        "solve_time" => solve_time(model),
        "primal_objective_value" => z_primal,
        "dual_objective_value" => z_dual,
    )
end

function extract_result(opf::OPFModel{<:AbstractFormulation})
    result = Dict{String,Any}()
    result["meta"] = extract_metadata(opf)
    result["primal"] = extract_primal(opf)
    result["dual"] = extract_dual(opf)
    return result
end
