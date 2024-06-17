# common distributions
include("glocal.jl")

abstract type AbstractOPFSampler end

function Random.rand(::AbstractRNG, ::AbstractOPFSampler)
    error("`rand` function not implemented for $(typeof(s)).")
end

struct SimpleOPFSampler{LS,RS}
    data::Dict
    load_sampler::LS
    reserve_sampler::RS
end

function SimpleOPFSampler(data::Dict, config::Dict)
    data = deepcopy(data)

    # Instantiate load sampler
    load_sampler = LoadScaler(data, config["load"])

    # Instantiate reserve sampler
    get!(config, "reserve", Dict())
    reserve_sampler = ReserveScaler(data, config["reserve"])

    return SimpleOPFSampler(data, load_sampler, reserve_sampler)
end

function Random.rand(rng::AbstractRNG, opf_sampler::SimpleOPFSampler)
    data = deepcopy(opf_sampler.data)
    rand!(rng, opf_sampler, data)
end

"""
    rand!(rng::AbstractRNG, s::AbstractOPFSampler, data::Dict)

Sample one new OPF instance and modify `data` in-place.

`data` must be a `Dict` in PowerModels format, representing the same network
    (i.e., same grid components with same indexing) as the one used to create `s`.
"""
function Random.rand!(rng::AbstractRNG, s::SimpleOPFSampler, data::Dict)
    pd, qd = rand(rng, s.load_sampler)
    _set_loads!(data, pd, qd)

    MRR, rmin, rmax, reserve_cost = rand(rng, s.reserve_sampler)
    _set_reserve!(data, MRR, rmin, rmax, reserve_cost)

    return data
end

function _set_loads!(data, pd, qd)
    L = length(data["load"])
    length(pd) == L || throw(DimensionMismatch())
    length(qd) == L || throw(DimensionMismatch())
    
    for i in 1:L
        ldat = data["load"]["$i"]
        ldat["pd"] = pd[i]
        ldat["qd"] = qd[i]
    end

    return nothing
end

function _set_reserve!(data, MRR, rmin, rmax, reserve_cost)
    G = length(data["gen"])
    length(rmin) == G || throw(DimensionMismatch())
    length(rmax) == G || throw(DimensionMismatch())

    for i in 1:length(data["gen"])
        gdat = data["gen"]["$i"]
        gdat["rmin"] = rmin[i]
        gdat["rmax"] = rmax[i]
        gdat["reserve_cost"] = reserve_cost[i, :]
    end

    data["minimum_reserve"] = MRR

    return nothing
end

include("load.jl")
include("reserve.jl")
