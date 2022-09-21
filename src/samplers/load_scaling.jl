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
    d::Uniform{Float64}

    # Reference active/reactive power demand
    pd_ref::Vector{Float64}
    qd_ref::Vector{Float64}

    # constructor
    function SimpleLoadScaling(pd::Vector{Float64}, qd::Vector{Float64}, l::Float64, u::Float64)
        L = length(pd)
        length(qd) == L || error("pd and qd must have same size")

        return new(Uniform(l, u), copy(pd), copy(qd))
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

function _sample_loads(rng, ls::SimpleLoadScaling)
    # Sample multiplicative noise and re-scale demand
    L = length(ls.pd_ref)
    α = rand(rng, ls.d, L)
    pd = ls.pd_ref .* α
    qd = ls.qd_ref .* α

    return pd, qd
end
