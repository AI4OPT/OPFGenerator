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
    σs2 = σs .^ 2
    d_η = MvLogNormal(-σs2 ./ 2, Diagonal(σs2))

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
