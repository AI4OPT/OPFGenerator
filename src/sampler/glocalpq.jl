"""
    GlocalPQ{G,L}

A glocal distribution with global/local factors `α::G` separate active/reactive factors `η₁::L` and `η₂::L`.

This distribution represents a random variable of the form `ϵ = α×η`, where
* `α` is a _scalar_ random variable, with distribution `d_α::G`
* `η₁` is a _vector_ random variable, with distribution `d_η::L`
* `η₂` is a _vector_ random variable, with distribution `d_η::L`
* `α`, `η₁`, and `η₂` are independent random variables
"""
struct GlocalPQ{G<:UvDist,L<:MvDist} <: MvDist
    d_α::G
    d_η::L
end

function GlocalPQ(d::Glocal)
    return GlocalPQ(d.d_α, d.d_η)
end

Distributions.length(d::GlocalPQ) = length(d.d_η)
Distributions.eltype(::GlocalPQ) = Float64

function Distributions._rand!(::AbstractRNG, ::GlocalPQ, ::AbstractArray)
    throw(DimensionMismatch(
        "Inconsistent argument dimensions: only vector and matrix-shaped `x` is supported."
    ))
end

function Distributions._rand!(rng::AbstractRNG, d::GlocalPQ, x::AbstractVector)
    length(x) == length(d) || throw(DimensionMismatch("Inconsistent argument dimensions."))

    x₁ = x
    x₂ = similar(x₁)

    α = rand(rng, d.d_α)
    η₁ = rand(rng, d.d_η)
    η₂ = rand(rng, d.d_η)

    x₁ .= α .* η₁
    x₂ .= α .* η₂

    return x₁, x₂
end

function Distributions._rand!(rng::AbstractRNG, d::GlocalPQ, x::AbstractMatrix)
    size(x, 1) == length(d) || throw(DimensionMismatch("Inconsistent argument dimensions."))
    n = size(x, 2)

    x₁ = x
    x₂ = similar(x₁)

    α = rand(rng, d.d_α, n)
    η₁ = rand(rng, d.d_η, n)
    η₂ = rand(rng, d.d_η, n)

    mul!(x₁, η₁, Diagonal(α))  # re-scales every column `η[:, j]` by `α[j]`
    mul!(x₂, η₂, Diagonal(α))

    return x₁, x₂
end
