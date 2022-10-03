"""
    SimpleLoadScaling

Load sampler that re-scales all loads by a multiplicative noise.

```math
\\begin{lmatrix}
    p^{d}_{i}\\
    q^{d}_{i}
\\end{lmatrix}
=
\\begin{lmatrix}
    \\alpha_{i} \\times \\bar{p}^{d}_{i}\\\\
    \\alpha_{i} \\times \\bar{q}^{d}_{i}
\\end{lmatrix}
```
where ``\\bar{p}^{d}, \\bar{q}^{d}`` are the reference active and reactive demands,
    and ``\\alpha_{i} \\sim U[l_{i}, u_{i}]``. The ``\\alpha_{i}``s are sampled i.i.d, i.e.,
    the noise is uncorrelated between loads.
"""
struct SimpleLoadScaling <: AbstractLoadSampler
    # Lower & upper bounds of multiplicative noise
    l::Vector{Float64}
    u::Vector{Float64}

    # Reference active/reactive power demand
    pd_ref::Vector{Float64}
    qd_ref::Vector{Float64}

    # Constructor
    function SimpleLoadScaling(l::Vector{Float64}, u::Vector{Float64}, pd::Vector{Float64}, qd::Vector{Float64})
        L = length(pd)
        length(qd) == L || error("pd and qd must have same size")
        length(l) == L || error("Lower bounds must be the same size as pd and qd")
        length(u) == L || error("Upper bounds must be the same size as pd and qd")

        all(l .<= u) || error("for all bounds (l, u), `l` must be smaller than or equal to `u`.")

        # Copy original data to avoid external changes afterwards
        return new(copy(l), copy(u), copy(pd), copy(qd))
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
    return SimpleLoadScaling(l .* ones(Float64, L), u .* ones(Float64, L), pd, qd)
end

function SimpleLoadScaling(data::Dict, l::Vector{Float64}, u::Vector{Float64})
    L = length(data["load"])
    pd = zeros(Float64, L)
    qd = zeros(Float64, L)
    for i in 1:L
        ldata = data["load"]["$i"]
        pd[i] = ldata["pd"]
        qd[i] = ldata["qd"]
    end
    return SimpleLoadScaling(l, u, pd, qd)
end

function _sample_loads(rng, ls::SimpleLoadScaling)
    L = length(ls.pd_ref)
    
    # Sample multiplicative noise
    d = Uniform(0.0, 1.0)
    α = rand(rng, d, L)
    α .*= (ls.u .- ls.l)
    α .+= ls.l  # this gives us α[i] ~ U(l[i], u[i])

    pd = ls.pd_ref .* α
    qd = ls.qd_ref .* α

    return pd, qd
end
