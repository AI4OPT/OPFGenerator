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

function Random.rand(rng::AbstractRNG, ls::LoadScaler{T}) where T <: Glocal
    ϵ = rand(rng, ls.d)  # sample multiplicative noise

    pd = ϵ .* ls.pd_ref
    qd = ϵ .* ls.qd_ref

    return pd, qd
end

function Random.rand(rng::AbstractRNG, ls::LoadScaler{T}) where T <: GlocalPQ
    ϵ₁, ϵ₂ = rand(rng, ls.d)

    pd = ϵ₁ .* ls.pd_ref
    qd = ϵ₂ .* ls.qd_ref

    return pd, qd
end

function LoadScaler(data::OPFData, options::Dict)
    L = data.L
    pd = data.pd
    qd = data.qd

    # Noise distribution
    # TODO: modularize this
    noise_type = get(options, "noise_type", "")
    if noise_type ∉ ["ScaledLogNormal", "ScaledUniform", "ScaledLogNormalPQ", "ScaledUniformPQ"]
        error("""Invalid noise type for load: $(noise_type). Supported values are:
        * \"ScaledLogNormal\"
        * \"ScaledUniform\"""")
    end
    
    # Grab noise parameters, and perform sanity checks
    l = get(options, "l", NaN)
    u = get(options, "u", NaN)
    if !(isa(l, Real) && isfinite(l))
        error("Missing or invalid input data: \"l\"")
    end
    if !(isa(u, Real) && isfinite(u))
        error("Missing or invalid input data: \"u\"")
    end
    if !(l <= u)
        error("Invalid global scaling parameters: [$l, $u]")
    end

    σs = zeros(Float64, L)
    σ = get(options, "sigma", NaN)
    if isa(σ, AbstractVector{<:Real})
        length(σ) == L || error("Invalid sigma length")
        σs .= σ
    elseif isa(σ, Real) && isfinite(σ)
        σs .= σ
    else
        error("Invalid input data: sigma. Must be a finite real number or a real-valued vector.")
    end

    # Noise distribution
    d = if noise_type == "ScaledLogNormal"
        ScaledLogNormal(l, u, σs)
    elseif noise_type == "ScaledUniform"
        ScaledUniform(l, u, σs)
    elseif noise_type == "ScaledLogNormalPQ"
        GlocalPQ(ScaledLogNormal(l, u, σs))
    elseif noise_type == "ScaledUniformPQ"
        GlocalPQ(ScaledUniform(l, u, σs))
    end

    return LoadScaler(d, pd, qd)
end
