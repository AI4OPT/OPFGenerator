# Contains a list of all supported OPF models
const SUPPORTED_OPF_MODELS = Type{<:AbstractPowerModel}[
    PowerModels.ACPPowerModel,
    PowerModels.DCPPowerModel,
    PowerModels.SOCWRPowerModel,
    PowerModels.SOCWRConicPowerModel,
]

# A name --> type dictionary
# Used for passing the OPF type as a string (e.g. through config file)
const OPF2TYPE = Dict{String,Type{<:Union{AbstractPowerModel,StandardFormOPFModel}}}(
    "ACPPowerModel" => PowerModels.ACPPowerModel,
    "DCPPowerModel" => PowerModels.DCPPowerModel,
    "SOCWRPowerModel" => PowerModels.SOCWRPowerModel,
    "SOCWRConicPowerModel" => PowerModels.SOCWRConicPowerModel,
    "StandardFormDCPPowerModel" => StandardFormDCPPowerModel,
)

include("acp.jl")    # ACCPPowerModel
include("dcp.jl")    # DCPPowerModel
include("socwr.jl")  # SOCWRPowerModel & SOCWRConicPowerModel

include("dcp_std.jl") # StandardFormOPFModel{DCPPowerModel}