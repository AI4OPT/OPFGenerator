using MathOptSymbolicAD

mutable struct OPFModel{OPF <: PM.AbstractPowerModel}
    data::Dict{String,Any}
    model::JuMP.GenericModel
end

include("acp.jl")      # ACPPowerModel
include("dcp.jl")      # DCPPowerModel
include("ed.jl")       # EconomicDispatch
include("socwr.jl")    # SOCWRPowerModel & SOCWRConicPowerModel
include("sdpwrm.jl")

# Contains a list of all supported OPF models
const SUPPORTED_OPF_MODELS = Type{<:AbstractPowerModel}[
    PowerModels.ACPPowerModel,
    PowerModels.DCPPowerModel,
    EconomicDispatch,
    PowerModels.SOCWRPowerModel,
    PowerModels.SOCWRConicPowerModel,
    PowerModels.SDPWRMPowerModel
]

# A name --> type dictionary
# Used for passing the OPF type as a string (e.g. through config file)
const OPF2TYPE = Dict{String,Type{<:AbstractPowerModel}}(
    "ACPPowerModel" => PowerModels.ACPPowerModel,
    "DCPPowerModel" => PowerModels.DCPPowerModel,
    "EconomicDispatch" => EconomicDispatch,
    "SOCWRPowerModel" => PowerModels.SOCWRPowerModel,
    "SOCWRConicPowerModel" => PowerModels.SOCWRConicPowerModel,
    "SDPWRMPowerModel" => PowerModels.SDPWRMPowerModel
)

function solve!(opf::OPFModel{<:AbstractPowerModel})
    optimize!(opf.model; _differentiation_backend = MathOptSymbolicAD.DefaultBackend())
end
