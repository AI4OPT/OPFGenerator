const UvDist = Distributions.ContinuousUnivariateDistribution
const MvDist = Distributions.ContinuousMultivariateDistribution

"""
    Glocal{G,L}

A glocal distribution with global/local factors `α::G` and  `η::L`.

This distribution represents a random variable of the form `ϵ = α×η`, where
* `α` is a _scalar_ random variable, with distribution `d_α::G`
* `η` is a _vector_ random variable, with distribution `d_η::L`
* `α` and `η` are independent random variables
"""
struct Glocal{G<:UvDist,L<:MvDist} <: MvDist
    d_α::G
    d_η::L
end

Distributions.length(d::Glocal) = length(d.d_η)
Distributions.eltype(::Glocal) = Float64

# custom samplers for Array-shaped data
function Distributions._rand!(::AbstractRNG, ::Glocal, ::AbstractArray)
    throw(DimensionMismatch(
        "Inconsistent argument dimensions: only vector and matrix-shaped `x` is supported."
    ))
end
function Distributions._rand!(rng::AbstractRNG, d::Glocal, x::AbstractVector)
    length(x) == length(d) || throw(DimensionMismatch("Inconsistent argument dimensions."))

    α = rand(rng, d.d_α)
    η = rand(rng, d.d_η)

    x .= α .* η

    return x
end
function Distributions._rand!(rng::AbstractRNG, d::Glocal, x::AbstractMatrix)
    size(x, 1) == length(d) || throw(DimensionMismatch("Inconsistent argument dimensions."))
    n = size(x, 2)

    α = rand(rng, d.d_α, n)
    η = rand(rng, d.d_η, n)

    mul!(x, η, Diagonal(α))  # re-scales every column `η[:, j]` by `α[j]`

    return x
end

"""
    ScaledLogNormal(l, u, σs)

Generate a [`Glocal`](@ref) distribution `ϵ = α×η` where `α~U[l,u]` and `ηᵢ ~ LogNormal(-σᵢ²/2, σᵢ)`.
"""
function ScaledLogNormal(l::Float64, u::Float64, σs::Vector{Float64})
    l <= u || error("Invalid bounds: l > u")
    all(σs .>= 0) || error("Invalid input: σs must be non-negative")

    # Uniform distribution
    d_α = Uniform(l, u)
    
    # Generate LogNormal distributions
    σs2 = σs .^ 2
    d_η = MvLogNormal(-σs2 ./ 2, Diagonal(σs2))

    return Glocal(d_α, d_η)
end

"""
    ScaledUniform(l, u, σs)

Generate a [`Glocal`](@ref) distribution `ϵ = α×η` where `α ~ U[l,u]` and `ηᵢ ~ U[1-σᵢ, 1+σᵢ]`.
"""
function ScaledUniform(l::Float64, u::Float64, σs::Vector{Float64})
    l <= u || error("Invalid bounds: l > u")
    all(σs .>= 0) || error("Invalid input: σs must be non-negative")

    # Uniform distribution
    d_α = Uniform(l, u)
    
    # Load-level uniform distributions
    d_η = Distributions.product_distribution([Uniform(1 - σ, 1 + σ) for σ in σs])

    return Glocal(d_α, d_η)
end
