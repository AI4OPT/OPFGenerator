module ACOPFGenerator

using Random
using StableRNGs
using Distributions
using Ipopt
using JuMP

using CodecBzip2
using CodecZlib
using JSON

using PowerModels
const PM = PowerModels

export load_json, save_json
export SimpleOPFSampler, SimpleLoadScaling, ScaleLogNorm
export sample
export sample_solve

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

"""
    sample_solve(sampler::AbstractOPFSampler, rng::StableRNG)

Adds noise to data in `sampler` based upon sampler Struct type and random seed `rng`.
"""
function sample_solve(rng::AbstractRNG, sampler::SimpleOPFSampler)
    agmtd_data = sample(rng, sampler)
    sol = solve_ac_opf(agmtd_data, Ipopt.Optimizer)
    meta = Dict("network" => sampler.data["name"],
            "seed" => rng,
            "augment_method" => typeof(sampler.load_sampler))
    return Dict("agmtd_data" => agmtd_data, "sol" => sol, "meta" => meta)
end

end # module
