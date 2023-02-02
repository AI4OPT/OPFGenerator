import Distributions: length, eltype, _rand!

struct LoadScaler{D} <: AbstractLoadSampler
    d::D
    pd_ref::Vector{Float64}
    qd_ref::Vector{Float64}
end

function _sample_loads(rng, ls::LoadScaler)
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

"""
    ScaledLogNormal

Generate correlated, multiplicative noise of the form `ϵ = α×η`, where
    ``α~U[l,u]`` is the same for all loads, and ``η~LogNormal(μ, σ)`` (uncorrelated).
The resulting multiplicative noise at load `i` has distribution `LogNormal(μᵢ+log α, σᵢ)`.

This distribution allows to control the total demand (via `α`),
    and the per-load noise level (via `σ`).
"""
struct ScaledLogNormal <: Distributions.ContinuousMultivariateDistribution
    d_α::Uniform{Float64}
    d_η::MvLogNormal
end

"""
    ScaledLogNormal(l, u, σs)

Generate a `ScaledLogNormal` distribution where `α~U[l,u]` and `ηᵢ ~ LogNormal(-σᵢ²/2, σᵢ)`.
"""
function ScaledLogNormal(l::Float64, u::Float64, σs::Vector{Float64})
    l <= u || error("Invalid bounds: l > u")

    # Uniform distribution
    d_α = Uniform(l, u)
    
    # Generate LogNormal distributions
    μs = - (σs .^ 2) ./ 2.0
    d_η = MvLogNormal(μs, σs)

    return ScaledLogNormal(d_α, d_η)
end

Distributions.length(d::ScaledLogNormal) = length(d.d_η)
Distributions.eltype(::ScaledLogNormal) = Float64

# custom samplers for Vector- and Matrix-shaped data
function Distributions._rand!(::AbstractRNG, ::ScaledLogNormal, ::AbstractArray)
    throw(DimensionMismatch(
        "Inconsistent argument dimensions: only vector and matrix-shaped `x` is supported."
    ))
end
function Distributions._rand!(rng::AbstractRNG, d::ScaledLogNormal, x::AbstractVector)
    length(x) == length(d) || throw(DimensionMismatch("Inconsistent argument dimensions."))

    α = rand(rng, d.d_α)
    η = rand(rng, d.d_η)

    x .= α .* η

    return x
end
function Distributions._rand!(rng::AbstractRNG, d::ScaledLogNormal, x::AbstractMatrix)
    size(x, 1) == length(d) || throw(DimensionMismatch("Inconsistent argument dimensions."))
    n = size(x, 2)

    α = rand(rng, d.d_α, n)
    η = rand(rng, d.d_η, n)

    for j in 1:n
        @views x[:, j] .= α[j] .* η[:, j]
    end

    return x
end
