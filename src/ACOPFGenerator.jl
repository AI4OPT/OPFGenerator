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

export load_json, save_json
export SimpleOPFSampler, SimpleLoadScaling, ScaleLogNorm
export sample

abstract type AbstractOPFSampler end

function sample(::AbstractRNG, ::AbstractOPFSampler)
    error("`sample` function not implemented for $(typeof(s)).")
end

abstract type AbstractLoadSampler end

include("utils.jl")

include("samplers/load_scaling.jl")
include("samplers/scale_log_norm.jl")

struct SimpleOPFSampler{LS}
    data::Dict
    load_sampler::LS
end

function sample(rng::AbstractRNG, opf_sampler::SimpleOPFSampler)
    data = deepcopy(opf_sampler.data)
    pd, qd = _sample_loads(rng, opf_sampler.load_sampler)
    _set_loads!(data, pd, qd)
    return data
end

include("acopf.jl")

end # module
