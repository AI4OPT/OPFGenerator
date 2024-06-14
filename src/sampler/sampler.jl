abstract type AbstractOPFSampler end

function Random.rand(::AbstractRNG, ::AbstractOPFSampler)
    error("`rand` function not implemented for $(typeof(s)).")
end

struct SimpleOPFSampler{LS}
    data::Dict
    load_sampler::LS
end

function SimpleOPFSampler(data::Dict, config::Dict)
    data = deepcopy(data)

    # Instantiate load sampler
    load_sampler = LoadScaler(data, config["load"])
    # Instantiate other stuff
    # TODO

    return SimpleOPFSampler(data, load_sampler)
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

include("load.jl")