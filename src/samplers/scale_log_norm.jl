"""
    ScaleLogNorm

Load sampler that re-scales loads, then scales by an individual log normal noise.

```math
\\begin{lmatrix}
    p^d_{i}\\
    q^d_{i}
\\end{lmatrix}
=
\\begin{lmatrix}
    \\alpha \\times \\eta_{i} \\times \\bar{p}^d_{i}\\
    \\alpha \\times \\eta_{i} \\times \\bar{q}^d_{i}
\\end{lmatrix}
```
where ``\\alpha is sampled from \\sim U[l, u]``, ``ηᵢ ~ LogNormal[-σᵢ²/2, σᵢ]``, 
    and ``\\bar{p}^d, \\bar{q}^d`` are the reference active and reactive demands.
"""
struct ScaleLogNorm <: AbstractLoadSampler
    # Upper and lower Uniform bounds
    d0::Uniform{Float64}

    # LogNormal Distributions
    ds::Vector{LogNormal}  # TODO: use MvLogNormal with Diagonal Covariance instead

    # Reference active/reactive power demand
    pd_ref::Vector{Float64}
    qd_ref::Vector{Float64}

    # constructor
    function ScaleLogNorm(pd::Vector{Float64}, qd::Vector{Float64}, l::Float64, u::Float64, σs::Vector{Float64})
        L = length(pd)
        length(qd) == L || error("pd and qd must have same size")
        l < u || error("l must be less than u")

        # Generate LogNormal distributions
        μs = - (σs .^ 2) ./ 2.0
        ds = LogNormal.(μs, σs)

        return new(Uniform(l, u), ds, copy(pd), copy(qd))
    end
end

function ScaleLogNorm(data::Dict, l::Float64, u::Float64, σ::Vector{Float64})
    L = length(data["load"])
    pd = zeros(Float64, L)
    qd = zeros(Float64, L)
    for i in 1:L
        ldata = data["load"]["$i"]
        pd[i] = ldata["pd"]
        qd[i] = ldata["qd"]
    end
    return ScaleLogNorm(pd, qd, l, u, σ)
end
ScaleLogNorm(data::Dict, l::Float64, u::Float64, σ::Float64) = ScaleLogNorm(data, l, u, fill(σ, length(data["load"])))

function _sample_loads(rng, ls::ScaleLogNorm)
    # Base load factor
    α = rand(rng, ls.d0)  # random scale value

    # Multiplicative noise
    η = rand.(rng, ls.ds)

    pd = α .* η .* ls.pd_ref
    qd = α .* η .* ls.qd_ref

    return pd, qd
end
