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
    data::Dict
end

function Random.rand(rng::AbstractRNG, rs::Nminus1StatusSampler)
    p = rand(rng)

    if p < 0.5
        # drop a branch 
        is_bridge = bridges(rs.data)
        non_bridge = [e for (e, b) in is_bridge if !b]
        
        e = rand(rng, non_bridge)
        e_index = rs.data["branch"][e]["index"]

        br_status = ones(Int, length(rs.data["branch"]))
        br_status[e_index] = 0

        gen_status = ones(Int, length(rs.data["gen"]))

        return br_status, gen_status

    else
        # drop a generator
        G = length(rs.data["gen"])
        g = rand(rng, 1:G)

        gen_status = ones(Int, G)
        gen_status[g] = 0

        br_status = ones(Int, length(rs.data["branch"]))

        return br_status, gen_status

    end
end


struct FullStatusSampler <: AbstractStatusSampler
    data::Dict
end

function Random.rand(rng::AbstractRNG, rs::FullStatusSampler)
    return ones(Int, length(rs.data["branch"])), ones(Int, length(rs.data["gen"]))
end

function StatusSampler(data::Dict, options::Dict)
    get(data, "basic_network", false) || error(
        """Invalid data: network data must be in basic format.
        Call `make_basic_network(data)` before calling this function"""
    )

    status_type = get(options, "type", "Full")
    if status_type == "Nminus1"
        return Nminus1StatusSampler(data)
    elseif status_type == "Full"
        return FullStatusSampler(data)
    else
        error("Invalid status sampler type: $status_type. Only 'Nminus1' and 'Full' are supported.")
    end
end
