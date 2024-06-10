module OPFGenerator

using Random
using StableRNGs
using Distributions

using PowerModels
const PM = PowerModels
using PGLib
using JuMP

import Random: rand, rand!

export load_json, save_json, save_h5
export SimpleOPFSampler, LoadScaler, ScaledLogNormal

include("utils.jl")
include("utils/json.jl")
include("utils/hdf5.jl")

include("bridges.jl")
include("process.jl")

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

include("opf/opf.jl")

end # module
