mutable struct OPFModel{OPF <: PM.AbstractPowerModel}
    data::Dict{String,Any}
    model::JuMP.Model
end

function build_opf(data, ::Type{<:PM.AbstractPowerModel}) end

include("acp.jl")    # ACCPPowerModel
include("dcp.jl")    # DCPPowerModel
include("socwr.jl")  # SOCWRPowerModel
