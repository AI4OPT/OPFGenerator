module ACOPFGenerator

using Random
using StableRNGs
using Distributions

using PowerModels
const PM = PowerModels
using JSON
using CodecBzip2

export SimpleOPFSampler, SimpleLoadScaling, sample, LNOPFSampler, ScaleLogNorm

abstract type AbstractOPFSampler end

function sample(::AbstractRNG, ::AbstractOPFSampler)
    error("`sample` function not implemented for $(typeof(s)).")
end

abstract type AbstractLoadSampler end

include("utils.jl")

include("samplers/load_scaling.jl")

include("samplers/scale_log_norm.jl")

struct SimpleOPFSampler
    data::Dict
    load_sampler::SimpleLoadScaling
end

struct LNOPFSampler
    data::Dict
    load_sampler::ScaleLogNorm
end

function sample(rng::AbstractRNG, opf_sampler::SimpleOPFSampler)
    data = deepcopy(opf_sampler.data)
    pd, qd = _sample_loads(rng, opf_sampler.load_sampler)
    _set_loads!(data, pd, qd)
    return data
end

function sample(rng::AbstractRNG, opf_sampler::LNOPFSampler)
    data = deepcopy(opf_sampler.data)
    pd, qd = _sample_loads(rng, opf_sampler.load_sampler)
    _set_loads!(data, pd, qd)
    return data
end

end # module
