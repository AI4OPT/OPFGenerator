using MathOptSymbolicAD

mutable struct OPFModel{OPF <: PM.AbstractPowerModel}
    data::Dict{String,Any}
    model::JuMP.Model
end

include("acp.jl")      # ACPPowerModel
include("dcp.jl")      # DCPPowerModel
include("dcp_ptdf.jl") # DCPPowerModel with PTDF
include("socwr.jl")    # SOCWRPowerModel & SOCWRConicPowerModel

# Contains a list of all supported OPF models
const SUPPORTED_OPF_MODELS = Type{<:AbstractPowerModel}[
    PowerModels.ACPPowerModel,
    PowerModels.DCPPowerModel,
    DCPPTDFPowerModel,
    PowerModels.SOCWRPowerModel,
    PowerModels.SOCWRConicPowerModel,
]

# A name --> type dictionary
# Used for passing the OPF type as a string (e.g. through config file)
const OPF2TYPE = Dict{String,Type{<:AbstractPowerModel}}(
    "ACPPowerModel" => PowerModels.ACPPowerModel,
    "DCPPowerModel" => PowerModels.DCPPowerModel,
    "DCPPTDFPowerModel" => DCPPTDFPowerModel,
    "SOCWRPowerModel" => PowerModels.SOCWRPowerModel,
    "SOCWRConicPowerModel" => PowerModels.SOCWRConicPowerModel,
)

function solve!(opf::OPFModel{<:AbstractPowerModel})
    optimize!(opf.model; _differentiation_backend = MathOptSymbolicAD.DefaultBackend())
end
