using MathOptSymbolicAD
using Clarabel
using Quadmath

function patch_if_clarabel!(model)
    if JuMP.solver_name(model) == "Clarabel"
        @warn "Removing VectorizeBridge from Clarabel (see MathOptInterface.jl#2452)"

        T = typeof(model).parameters[1]

        JuMP.MOI.Bridges.remove_bridge(
           JuMP.backend(model).optimizer,
           JuMP.MOI.Bridges.Variable.VectorizeBridge{T},
       )
    end
end

mutable struct OPFModel{OPF <: PM.AbstractPowerModel}
    data::Dict{String,Any}
    model::JuMP.GenericModel
end

include("acp.jl")      # ACPPowerModel
include("dcp.jl")      # DCPPowerModel
include("socwr.jl")    # SOCWRPowerModel & SOCWRConicPowerModel

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

function solve!(opf::OPFModel{<:AbstractPowerModel})
    optimize!(opf.model; _differentiation_backend = MathOptSymbolicAD.DefaultBackend())
end
