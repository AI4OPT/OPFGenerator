import Distributions: length, eltype, _rand!

abstract type AbstractReserveSampler end

"""
    E2ELRReserveScaler{U,F}

Samples reserve requirements following the procedure below:

1. Sample a minimum reserve requirement `MRR` from a uniform distribution `U(lb, ub)`.
2. Compute the upper bound of reserve requirements for each generator as `rmax = α * (pmax - pmin)`.
3. Fix the lower bound of reserve requirement per generator to zero.

The parameter `α` is a scaling factor that determines the total reserve requirement
    as a multiple of the maximum power range of all generators.

"""
struct E2ELRReserveScaler <: AbstractReserveSampler
    mrr_dist::Uniform
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

# """
#     UCjlReserveScaler

# Samples reserve requirements following the procedure below:

# 1. Compute minimum reserve requirement `MRR` by scaling the reference power demand
#     with the minimum reserve fraction.

# 2. Compute the upper bound of reserve requirements for each generator by scaling
#     the active generation upper bound with the maximum reserve fraction.

# NOTE: This reserve scaler itself does not use any random variables to compute
#     reserve requirements. However, it depends on the active load and generation
#     upper bound, which themselves may vary for each sample.

# """
# struct UCjlReserveScaler <: AbstractReserveSampler
#     min_frac::Float64
#     max_frac::Float64
#     pd_ref::Vector{Float64}
#     pg_max::Vector{Float64}
# end

# function Random.rand(rng::AbstractRNG, rs::UCjlReserveScaler)
#     # generate MRR
#     MRR = rs.min_frac * sum(rs.pd_ref)

#     rmax = rs.max_frac .* rs.pg_max
#     rmin = zeros(Float64, length(rmax))

#     return MRR, rmin, rmax
# end

function ReserveScaler(data::Dict, options::Dict)
    get(data, "basic_network", false) || error(
        """Invalid data: network data must be in basic format.
        Call `make_basic_network(data)` before calling this function"""
    )

    reserve_type = get(options, "type", "")
    if reserve_type == "E2ELR"
        l = options["l"]
        u = options["u"]
        factor = get(options, "factor", 5.0)
        mrr_dist = Uniform(l, u)

        pg_min = [gen["pmin"] for (id, gen) in data["gen"]]
        pg_max = [gen["pmax"] for (id, gen) in data["gen"]]

        return E2ELRReserveScaler(mrr_dist, factor, pg_min, pg_max)

    # elseif reserve_type == "UCjl"
    #     min_frac = get(options, "min_frac", 0.1)
    #     max_frac = get(options, "max_frac", 0.1)

    #     pd_ref = [data["load"]["$i"]["pd"] for i in 1:length(data["load"])]
    #     pg_max = [gen["pmax"] for (id, gen) in data["gen"]]

    #     @warn "When using the UCjl reserve scaler, the minimum reserve requirement (MRR)
    #          depends on the total active load, which typically varies per sample. Proceed with caution."

    #     return UCjlReserveScaler(min_frac, max_frac, pd_ref, pg_max)
    else
        error("Invalid noise type: $(reserve_type).\nOnly \"E2ELR\" and \"UCjl\" are supported.")
    end
end
