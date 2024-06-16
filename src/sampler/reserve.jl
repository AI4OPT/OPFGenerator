import Distributions: length, eltype, _rand!

abstract type AbstractReserveSampler end

"""
    E2ELRReserveScaler

Samples reserve requirements following the procedure below:

1. Sample a minimum reserve requirement `MRR` from a uniform distribution `U(lb, ub)` (`mrr_dist`).
2. Compute the upper bound of reserve requirements for each generator as `rmax = α * (pmax - pmin)`.
3. Fix the lower bound of reserve requirement per generator to zero.

The parameter `α` is a scaling factor that determines each generator's maximum reserve
    capacity. It is the `factor` parameter times the ratio of the largest generator's capacity
    to the sum of all generators' dispatchable capacity.

"""
struct E2ELRReserveScaler <: AbstractReserveSampler
    mrr_dist::Uniform{Float64}
    factor::Float64
    pg_min::Vector{Float64}
    pg_max::Vector{Float64}
end

function Random.rand(rng::AbstractRNG, rs::E2ELRReserveScaler)
    # generate MRR
    pmax = maximum(rs.pg_max)
    MRR = rand(rng, rs.mrr_dist) * pmax

    # generate reserve requirements
    pg_ranges = max.(0.0, rs.pg_max .- rs.pg_min)
    α = rs.factor * pmax / sum(pg_ranges)
    rmax = α .* pg_ranges
    rmin = zeros(Float64, length(rmax))
    return MRR, rmin, rmax
end


struct NullReserveScaler <: AbstractReserveSampler
    n_gen::Int
end

function Random.rand(rng::AbstractRNG, rs::NullReserveScaler)
    return 0.0, zeros(Float64, rs.n_gen), zeros(Float64, rs.n_gen)
end

function ReserveScaler(data::Dict, options::Dict)
    get(data, "basic_network", false) || error(
        """Invalid data: network data must be in basic format.
        Call `make_basic_network(data)` before calling this function"""
    )

    reserve_type = get(options, "type", "")
    if reserve_type == "E2ELR"
        @assert haskey(options, "l") "Missing lower bound `l` in config for E2ELR reserve sampler"
        @assert haskey(options, "u") "Missing upper bound `u` in config for E2ELR reserve sampler"
        @assert haskey(options, "factor") "Missing scaling factor `factor` in config for E2ELR reserve sampler"

        l = options["l"]
        u = options["u"]
        factor = options["factor"]
        mrr_dist = Uniform(l, u)

        pg_min = [gen["pmin"] for (id, gen) in data["gen"]]
        pg_max = [gen["pmax"] for (id, gen) in data["gen"]]

        return E2ELRReserveScaler(mrr_dist, factor, pg_min, pg_max)
    elseif reserve_type == ""
        n_gen = length(data["gen"])
        return NullReserveScaler(n_gen)
    else
        error("Invalid noise type: $(reserve_type).\nOnly \"E2ELR\", or no reserves, is supported.")
    end
end
