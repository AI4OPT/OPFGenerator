"""
    SimpleLoadScaling

Load sampler that re-scales all loads by a multiplicative noise.

```math
\\begin{lmatrix}
    p^{d}_{i}\\
    q_{d}_{i}
\\end{lmatrix}
=
\\begin{lmatrix}
    \\alpha_{i} \\times \\bar{p}^{d}_{i}\\
    \\alpha_{i} \\times \\bar{q}_{d}_{i}
\\end{lmatrix}
```
where ``\\alpha_{i}``s are sampled i.i.d from \\sim U[l, u]``, and ``\\bar{p}^{d}, \\bar{q}^{d}``
    are the reference active and reactive demands.
"""
struct SimpleLoadScaling <: AbstractLoadSampler
    # d::Uniform{Float64} # Q about generalizing this feature, as well as having r take in only [0,1]
    r::Matrix{Float64}

    # Reference active/reactive power demand
    pd_ref::Vector{Float64}
    qd_ref::Vector{Float64}

    # constructor
    function SimpleLoadScaling(pd::Vector{Float64}, qd::Vector{Float64}, l::Float64, u::Float64)
        L = length(pd)
        length(qd) == L || error("pd and qd must have same size")

        # return new(Uniform(l, u), copy(pd), copy(qd))
        return new(hcat(fill(l, L), fill(u, L))', copy(pd), copy(qd))
    end
    
    function SimpleLoadScaling(pd::Vector{Float64}, qd::Vector{Float64}, a::Vector{Tuple{Float64, Float64}})
        L = length(pd)
        length(qd) == L || error("pd and qd must have same size")

        a_mtx = hcat(collect.(a)...)
        length(a_mtx[1, :]) == L || error("'a' must be the same size as pd and qd")

        all(a_mtx[1, :] .<= a_mtx[2, :]) || error("for all tuples (l, u), l <= u.")

        return new(a_mtx, copy(pd), copy(qd))
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

function SimpleLoadScaling(data::Dict, a::Vector{Tuple{Float64, Float64}})
    L = length(data["load"])
    pd = zeros(Float64, L)
    qd = zeros(Float64, L)
    for i in 1:L 
        ldata = data["load"]["$i"]
        pd[i] = ldata["pd"]
        qd[i] = ldata["qd"]
    end
    return SimpleLoadScaling(pd, qd, a)
end

function _sample_loads(rng, ls::SimpleLoadScaling)
    # Sample multiplicative noise and re-scale demand
    l = ls.r[1, :]
    u = ls.r[2, :]
    L = length(ls.pd_ref)
    # scaling vector
    α = Vector{Float64}(undef, L)
    # distribution vector
    # d_pd = Vector{Uniform{Float64}}(undef, L)
    # d_qd = Vector{Uniform{Float64}}(undef, L)
      
    for i in 1:L
        dist = Uniform(l[i], u[i])
        α[i] = rand(rng, dist)
        # d_pd[i] = dist*ls.pd_ref[i]
        # d_qd[i] = dist*ls.qd_ref[i]
    end
    pd = ls.pd_ref .* α
    qd = ls.qd_ref .* α

    return pd, qd#, d_pd, d_qd
end