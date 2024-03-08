using MathOptSymbolicAD
using Clarabel
using Quadmath

function patch_if_clarabel!(model, optimizer)
    if (
        (optimizer == Clarabel.Optimizer) ||
        (
            (optimizer isa MOI.OptimizerWithAttributes) && (
                (optimizer.optimizer_constructor == Clarabel.Optimizer) ||
                (optimizer.optimizer_constructor <: Clarabel.Optimizer))
        )
    )
        @warn "Removing VectorizeBridge from Clarabel (see MathOptInterface.jl#2452)"

        T = if (
            (optimizer == Clarabel.Optimizer{Float128}) ||
            (
                (optimizer isa MOI.OptimizerWithAttributes) && (
                    (optimizer.optimizer_constructor == Clarabel.Optimizer{Float128}) ||
                    (optimizer.optimizer_constructor <: Clarabel.Optimizer{Float128}))
            )
        )
            Float128
        else    
            Float64
        end

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
