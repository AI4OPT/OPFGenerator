mutable struct OPFModel{OPF <: PM.AbstractPowerModel}
    data::Dict{String,Any}
    model::JuMP.Model
end

function build_opf(data, ::Type{<:PM.AbstractPowerModel}) end

include("acopf.jl")  # ACCPPowerModel
include("dcopf.jl")  # DCPPowerModel
include("socwr.jl")  # SOCWRPowerModel
