module ACOPFGenerator

using Random
using StableRNGs
using Distributions

using CodecBzip2
using CodecZlib
using JSON

using PowerModels
const PM = PowerModels
using JuMP

import Random: rand

export load_json, save_json
export SimpleOPFSampler, LoadScaler, ScaledLogNormal

include("utils.jl")
include("bridges.jl")

# Load samplers
include("load/load.jl")

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
    # Sample and update loads
    pd, qd = rand(rng, opf_sampler.load_sampler)
    _set_loads!(data, pd, qd)
    return data
end

include("acopf.jl")

end # module
