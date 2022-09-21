module ACOPFGenerator

using Random
using StableRNGs
using Distributions

using PowerModels
const PM = PowerModels
using JSON
using CodecBzip2

export SimpleOPFSampler, SimpleLoadScaling, sample

abstract type AbstractOPFSampler end

function sample(::AbstractRNG, ::AbstractOPFSampler)
    error("`sample` function not implemented for $(typeof(s)).")
end

abstract type AbstractLoadSampler end

include("utils.jl")

include("samplers/load_scaling.jl")

struct SimpleOPFSampler
    data::Dict        # Original data dictionary
    load_sampler::SimpleLoadScaling  # Sampler for load 
end

function sample(rng::AbstractRNG, opf_sampler::SimpleOPFSampler)
    data = deepcopy(opf_sampler.data)
    pd, qd = _sample_loads(rng, opf_sampler.load_sampler)
    _set_loads!(data, pd, qd)
    return data
end

end # module
