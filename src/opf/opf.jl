using MathOptSymbolicAD

# Contains a list of all supported OPF models
const SUPPORTED_OPF_MODELS = Type{<:AbstractPowerModel}[
    PowerModels.ACPPowerModel,
    PowerModels.DCPPowerModel,
    PowerModels.SOCWRPowerModel,
    PowerModels.SOCWRConicPowerModel,
]

# A name --> type dictionary
# Used for passing the OPF type as a string (e.g. through config file)
const OPF2TYPE = Dict{String,Type{<:AbstractPowerModel}}(
    "ACPPowerModel" => PowerModels.ACPPowerModel,
    "DCPPowerModel" => PowerModels.DCPPowerModel,
    "SOCWRPowerModel" => PowerModels.SOCWRPowerModel,
    "SOCWRConicPowerModel" => PowerModels.SOCWRConicPowerModel,
)

mutable struct OPFModel{OPF <: PM.AbstractPowerModel}
    data::Dict{String,Any}
    model::JuMP.Model
end

function solve!(opf::OPFModel{<:AbstractPowerModel})
    optimize!(opf.model; _differentiation_backend = MathOptSymbolicAD.DefaultBackend())
end

include("acp.jl")    # ACCPPowerModel
include("dcp.jl")    # DCPPowerModel
include("socwr.jl")  # SOCWRPowerModel & SOCWRConicPowerModel
