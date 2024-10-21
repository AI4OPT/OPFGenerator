import Distributions: length, eltype, _rand!

abstract type AbstractStatusSampler end

"""
    Nminus1StatusSampler

Samples a branch or generator to be inactive following the procedure below:

1. With probability 1/2, decide if a branch or generator is inactive in this instance.

If a branch is inactive:
2. Identify the set of branches which are not bridges. Removing a bridge results in a disconnected network.
3. Sample a branch to be inactive from the set of non-bridge branches.
4. Construct the status vector to be all ones except for zero at the sampled branch.

If a generator is inactive:
2. Sample a generator to be inactive.
3. Construct the status vector to be all ones except for zero at the sampled generator.
"""
struct Nminus1StatusSampler <: AbstractStatusSampler
    E::Int
    G::Int
    non_bridges::Vector{Int}
end

function Random.rand(rng::AbstractRNG, rs::Nminus1StatusSampler)
    p = rand(rng)

    br_status = ones(Int, rs.E)
    gen_status = ones(Int, rs.G)

    if p < 0.5
        # drop a branch 
        e = rand(rng, rs.non_bridges)
        br_status[e] = 0
    else
        # drop a generator
        g = rand(rng, 1:rs.G)
        gen_status[g] = 0
    end

    return br_status, gen_status
end


struct FullStatusSampler <: AbstractStatusSampler
    E::Int
    G::Int
end

function Random.rand(rng::AbstractRNG, rs::FullStatusSampler)
    return ones(Int, rs.E), ones(Int, rs.G)
end

function StatusSampler(data::OPFData, options::Dict)
    status_type = lowercase(get(options, "type", ""))
    if status_type == "nminus1"
        is_bridge = bridges(data)
        non_bridges = [e for (e, b) in enumerate(is_bridge) if !b]
        return Nminus1StatusSampler(data.E, data.G, non_bridges)
    elseif status_type == "full" || status_type == ""
        return FullStatusSampler(data.E, data.G)
    else
        error("Invalid status sampler type: $status_type. Only 'Nminus1' and 'Full' are supported.")
    end
end