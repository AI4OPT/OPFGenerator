import Distributions: length, eltype, _rand!

abstract type AbstractLoadSampler end

"""
    LoadScaler{D}

Scales loads with multiplicative noise sample from `d::D`.

Scaling retains power factors, i.e., each load's active 
    and reactive demand is scaled by the same number.
"""
struct LoadScaler{D} <: AbstractLoadSampler
    d::D
    pd_ref::Vector{Float64}
    qd_ref::Vector{Float64}
end

function Random.rand(rng::AbstractRNG, ls::LoadScaler)
    ϵ = rand(rng, ls.d)  # sample multiplicative noise

    pd = ϵ .* ls.pd_ref
    qd = ϵ .* ls.qd_ref

    return pd, qd
end

function LoadScaler(data::Dict, options::Dict)
    get(data, "basic_network", false) || error(
        """Invalid data: network data must be in basic format.
        Call `make_basic_network(data)` before calling this function"""
    )
    L = length(data["load"])
    pd = [data["load"]["$i"]["pd"] for i in 1:L]
    qd = [data["load"]["$i"]["qd"] for i in 1:L]

    # Noise distribution
    # TODO: modularize this
    noise_type = get(options, "noise_type", "")
    noise_type == "ScaledLogNormal" || error(
        "Invalid noise type: $(noise_type).\nOnly \"ScaledLogNormal\" is supported."
    )
    l = options["l"]
    u = options["u"]
    σ = options["sigma"]
    σs = zeros(Float64, L)
    if isa(σ, Vector)
        length(σ) == L || error("Invalid sigma length")
        σs .= σ
    elseif isa(σ, Real)
        σs .= σ
    else
        error("Invalid input data: sigma has type $(typeof(σ)) (must be a number or a float)")
    end
    d = ScaledLogNormal(l, u, σs)
    
    return LoadScaler(d, pd, qd)
end
