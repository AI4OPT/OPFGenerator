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
where ``\\alpha is sampled from \\sim U[l, u]``, \\eta_{i}``s 
    are sampled i.i.d from \\sim LogNormal[1, \\epsilon_{i}], 
    and ``\\bar{p}^d, \\bar{q}^d`` are the reference active and reactive demands.
"""
struct ScaleLogNorm <: AbstractLoadSampler
    # Upper and lower Uniform bounds
    α::Tuple{Float64, Float64}

    # Standard deviations of i.i.d. LogNormal distributions
    # size (1, L) where L = Float64 of loads
    σ::Vector{Float64}

    # Reference active/reactive power demand
    pd_ref::Vector{Float64}
    qd_ref::Vector{Float64}

    # constructor
    function ScaleLogNorm(pd::Vector{Float64}, qd::Vector{Float64}, l::Float64, u::Float64, σ::Vector{Float64}) 
        L = length(pd)
        length(qd) == L || error("pd and qd must have same size")
        l < u || error("l must be less than u")

        return new((l, u), copy(σ), copy(pd), copy(qd))
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

function _sample_loads(rng, ls::ScaleLogNorm)
    # Sample multiplicative noise and re-scale demand
    scale = rand(rng, Uniform(ls.α[1], ls.α[2])) # random scale value
    L = length(ls.pd_ref)

    # lognormal sample vector
    logn = Vector{Float64}(undef, L)

    for i in 1:L
        dist = LogNormal(1, ls.σ[i])
        logn[i] = rand(rng, dist)
    end

    pd = ls.pd_ref .* logn .* scale
    qd = ls.qd_ref .* logn .* scale
    return pd, qd
end