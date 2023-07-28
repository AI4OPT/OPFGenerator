# Contains a list of all supported OPF models
# This array is initialized as emtpy
const SUPPORTED_OPF_MODELS = Type{<:AbstractPowerModel}[
    PowerModels.ACPPowerModel,
    PowerModels.DCPPowerModel,
    PowerModels.SOCWRPowerModel,
    PowerModels.SOCWRConicPowerModel,
]

# A name --> type dictionary
# Used for passing the OPF type as a string (e.g. through config file)
const OPF2TYPE = Dict{String,Type{<:AbstractPowerModel}}(
    string(OPF) => OPF
    for OPF in SUPPORTED_OPF_MODELS
)

mutable struct OPFModel{OPF <: PM.AbstractPowerModel}
    data::Dict{String,Any}
    model::JuMP.Model
end

function build_opf(data, ::Type{<:PM.AbstractPowerModel}) end

include("acp.jl")    # ACCPPowerModel
include("dcp.jl")    # DCPPowerModel
include("socwr.jl")  # SOCWRPowerModel & SOCWRConicPowerModel
