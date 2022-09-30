"""
    SimpleLoadScaling

Load sampler that re-scales all loads by a multiplicative noise.

```math
\\begin{lmatrix}
    pd_{i}\\
    qd_{i}
\\end{lmatrix}
=
\\begin{lmatrix}
    \\alpha_{i} \\times \\bar{pd}_{i}\\
    \\alpha_{i} \\times \\bar{qd}_{i}
\\end{lmatrix}
```
where ``\\alpha_{i}``s are sampled i.i.d from \\sim U[l, u]``, and ``\\bar{pd}, \\bar{qd}``
    are the reference active and reactive demands.
"""
struct SimpleLoadScaling <: AbstractLoadSampler
    # Upper and lower Uniform bounds, size (2, L) where L = number of loads
    r::Matrix{Float64}

    # Reference active/reactive power demand
    pd_ref::Vector{Float64}
    qd_ref::Vector{Float64}

    # constructor
    function SimpleLoadScaling(pd::Vector{Float64}, qd::Vector{Float64}, l::Float64, u::Float64)
        L = length(pd)
        length(qd) == L || error("pd and qd must have same size")

        return new(hcat(fill(l, L), fill(u, L))', copy(pd), copy(qd))
    end
    
    function SimpleLoadScaling(pd::Vector{Float64}, qd::Vector{Float64}, r::Matrix{Float64})
        L = length(pd)
        length(qd) == L || error("pd and qd must have same size")
        length(r[1, :]) == L || error("'r' must be the same size as pd and qd")
        all(r[1, :] .< r[2, :]) || error("for all bounds (l, u), l must be less than u.")

        return new(r, copy(pd), copy(qd))
    end
    
end

function SimpleLoadScaling(data::Dict, l::Float64, u::Float64)
    L = length(data["load"])
    pd = zeros(Float64, L)
    qd = zeros(Float64, L)
    for i in 1:L
        ldata = data["load"]["$i"]
        pd[i] = ldata["pd"]
        qd[i] = ldata["qd"]
    end
    return SimpleLoadScaling(pd, qd, l, u)
end

function SimpleLoadScaling(data::Dict, r::Matrix{Float64})
    L = length(data["load"])
    pd = zeros(Float64, L)
    qd = zeros(Float64, L)
    for i in 1:L 
        ldata = data["load"]["$i"]
        pd[i] = ldata["pd"]
        qd[i] = ldata["qd"]
    end
    return SimpleLoadScaling(pd, qd, r)
end

function _sample_loads(rng, ls::SimpleLoadScaling)
    # Sample multiplicative noise and re-scale demand
    l = ls.r[1, :]
    u = ls.r[2, :]
    L = length(ls.pd_ref)
    # scaling vector
    α = Vector{Float64}(undef, L)

    for i in 1:L
        dist = Uniform(l[i], u[i])
        α[i] = rand(rng, dist)
    end

    pd = ls.pd_ref .* α
    qd = ls.qd_ref .* α

    return pd, qd
end