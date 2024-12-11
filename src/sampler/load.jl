import Distributions: length, eltype, _rand!

abstract type AbstractLoadSampler end

"""
    LoadScaler{D}

Scales loads with multiplicative noise sampled from `d::D`.

The distribution `d::D` is a `2*L`-dimensional distribution.
    The sample active/reactive demand for load ``i`` is denoted by `\\tilde{p}_{i}, \\tilde{q}_{i}`` and has the form
    ``\\tilde{p}_{i} = \\epsilon_{i} \\bar{p}_{i}``,
    ``\\tilde{q}_{i} = \\epsilon_{i+L} \\bar{q}_{i}``,
    where
    * ``\\bar{p}_{i}, \\bar{q}_{i}`` are the reference active/reactive demand for load ``i``
    * ``\\epsilon \\in \\mathbb{R}^{2L}`` is multiplicative noise sampled from distribution `d::D`.
"""
struct LoadScaler{D} <: AbstractLoadSampler
    d::D    # distribution of multiplicative noise
    pd_ref::Vector{Float64}  # reference active load
    qd_ref::Vector{Float64}  # reference reactive load
end

function Random.rand(rng::AbstractRNG, ls::LoadScaler)
    L = length(ls.pd_ref)
    ϵ = rand(rng, ls.d)  # sample multiplicative noise

    pd = ϵ[1:L] .* ls.pd_ref
    qd = ϵ[(L+1):(2*L)] .* ls.qd_ref

    return pd, qd
end

function LoadScaler(data::OPFData, options::Dict)
    L = data.L
    pd = data.pd
    qd = data.qd

    # Noise distribution
    # TODO: modularize this
    noise_type = get(options, "noise_type", "")
    if noise_type ∉ ["ScaledLogNormal", "ScaledUniform"]
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

    σs = zeros(Float64, 2*L)
    σ = get(options, "sigma", NaN)
    if isa(σ, AbstractVector{<:Real})
        length(σ) == L || error("Invalid sigma length")
        # We use the same noise level for active/reactive load
        σs[1:L] .= σ
        σs[(L+1):(2*L)] .= σ
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
    end

    return LoadScaler(d, pd, qd)
end
