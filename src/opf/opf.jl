using MathOptSymbolicAD
using SparseArrays

abstract type Formulation end

mutable struct OPFModel{OPF <: Formulation}
    data::Dict{String,Any}
    model::JuMP.GenericModel
end

struct OPFData
    case::String
    base_mva::Float64

    N::Int  # number of buses
    E::Int  # number of branches
    G::Int  # number of generators
    L::Int  # number of loads
    
    A::SparseMatrixCSC{Float64,Int}  # from/to branch incidence matrix
    Ag::SparseMatrixCSC{Float64,Int}  # generator incidence matrix

    # Bus data
    vnom::Vector{Float64}
    vmin::Vector{Float64}
    vmax::Vector{Float64}
    gs::Vector{Float64}  # shunt
    bs::Vector{Float64}  # shunt
    pd::Vector{Float64}  # nodal active power load
    qd::Vector{Float64}  # nodal reactive power load
    bus_arcs_fr::Vector{Vector{Int}}  # indices of branches exiting the bus
    bus_arcs_to::Vector{Vector{Int}}  # indices of branches entering the bus
    bus_gens::Vector{Vector{Int}}  # indices of generators at the bus
    ref_bus::Int  # index of slack bus

    # Generator data
    pgmin::Vector{Float64}
    pgmax::Vector{Float64}
    qgmin::Vector{Float64}
    qgmax::Vector{Float64}
    c0::Vector{Float64}  # constant cost
    c1::Vector{Float64}  # linear cost
    c2::Vector{Float64}  # quadratic cost
    gen_status::Vector{Bool}  # generator status

    # Branch data
    bus_fr::Vector{Int}  # from bus
    bus_to::Vector{Int}  # to bus
    g::Vector{Float64}
    b::Vector{Float64}
    gff::Vector{Float64}
    gft::Vector{Float64}
    gtf::Vector{Float64}
    gtt::Vector{Float64}
    bff::Vector{Float64}
    bft::Vector{Float64}
    btf::Vector{Float64}
    btt::Vector{Float64}
    smax::Vector{Float64}  # max apparent power flow
    dvamin::Vector{Float64}  # angle difference min
    dvamax::Vector{Float64}  # angle difference max
    branch_status::Vector{Bool}  # branch status
end

"""
    OPFData(network::Dict{String,Any})

Convert a PowerModels data dictionary to `OPFData` structure.

The PowerModels data dictionary must be in basic format.
"""
function OPFData(network::Dict{String,Any})
    @assert network["basic_network"] "Network data is not in basic format."
    @assert network["per_unit"] == true "Network data is not per-unit scaled."

    N = length(network["bus"])
    E = length(network["branch"])
    G = length(network["gen"])
    L = length(network["load"])

    # Bus data
    vnom = [network["bus"]["$i"]["base_kv"] for i in 1:N]
    vmin = [network["bus"]["$i"]["vmin"] for i in 1:N]
    vmax = [network["bus"]["$i"]["vmax"] for i in 1:N]

    # Aggregate shunts at the bus level
    gs = zeros(Float64, N)
    bs = zeros(Float64, N)
    for (_, shunt) in network["shunt"]
        shunt["status"] == 1 || continue  # skip inactive buses
        i = shunt["shunt_bus"]
        gs[i] += shunt["gs"]
        bs[i] += shunt["bs"]
    end

    # Aggregate loads at the bus level
    pd = zeros(Float64, N)
    qd = zeros(Float64, N)
    for (_, load) in network["load"]
        load["status"] == 1 || continue  # skip inactive loads
        i = load["load_bus"]
        pd[i] += load["pd"]
        qd[i] += load["qd"]
    end

    # Reference bus
    ref_buses = [i for i in 1:N if network["bus"]["$i"]["bus_type"] == 3]
    @assert length(ref_buses) == 1 "There must be exactly one reference bus"
    ref_bus = ref_buses[1]

    # Generator data
    pgmin = zeros(Float64, G)
    pgmax = zeros(Float64, G)
    qgmin = zeros(Float64, G)
    qgmax = zeros(Float64, G)
    c0 = zeros(Float64, G)
    c1 = zeros(Float64, G)
    c2 = zeros(Float64, G)
    gen_status = zeros(Bool, G)
    bus_gens = [Int[] for _ in 1:N]
    Ag_i = zeros(Int, G)
    Ag_j = zeros(Int, G)
    Ag_v = zeros(Float64, G)
    for g in 1:G
        gen = network["gen"]["$g"]

        i = gen["gen_bus"]
        push!(bus_gens[i], g)

        pgmin[g] = gen["pmin"]
        pgmax[g] = gen["pmax"]
        qgmin[g] = gen["qmin"]
        qgmax[g] = gen["qmax"]
        # ⚠️ cost data assumes quadratic cost everywhere
        c0[g] = gen["cost"][3]
        c1[g] = gen["cost"][2]
        c2[g] = gen["cost"][1]

        gen_status[g] = gen["gen_status"] == 1

        # Generator incidence matrix
        Ag_i[g] = network["gen"]["$g"]["gen_bus"]
        Ag_j[g] = g
        Ag_v[g] = network["gen"]["$g"]["gen_status"]
    end
    # sort everything again
    sort!.(bus_gens)

    Ag = sparse(Ag_i, Ag_j, Ag_v, N, G)

    # Branch data
    bus_fr = zeros(Int, E)
    bus_to = zeros(Int, E)
    branch_g = zeros(Float64, E)
    branch_b = zeros(Float64, E)
    gff = zeros(Float64, E)
    gft = zeros(Float64, E)
    gtf = zeros(Float64, E)
    gtt = zeros(Float64, E)
    bff = zeros(Float64, E)
    bft = zeros(Float64, E)
    btf = zeros(Float64, E)
    btt = zeros(Float64, E)
    smax = zeros(Float64, E)
    dvamin = zeros(Float64, E)
    dvamax = zeros(Float64, E)
    branch_status = zeros(Bool, E)
    bus_arcs_fr = [Int[] for _ in 1:N]
    bus_arcs_to = [Int[] for _ in 1:N]
    A_i = zeros(Int, 2*E)
    A_j = zeros(Int, 2*E)
    A_v = zeros(Float64, 2*E)
    for e in 1:E
        branch = network["branch"]["$e"]
        i::Int = branch["f_bus"]
        j::Int = branch["t_bus"]

        bus_fr[e] = i
        bus_to[e] = j
        push!(bus_arcs_fr[i], e)
        push!(bus_arcs_to[j], e)

        z::ComplexF64 = branch["br_r"] + im * branch["br_x"]
        y = inv(z)  # compute branch admittance
        isfinite(y) || error("Branch $e has zero impedance")
        g, b = real(y), imag(y)

        branch_g[e] = g
        branch_b[e] = b

        # Compute tap ratio
        τ::Float64 = get(branch, "tap", 1.0)
        ϕ::Float64 = get(branch, "shift", 0.0)
        tr = τ * cos(ϕ)
        ti = τ * sin(ϕ)
        ttm = abs2(τ)

        g_fr::Float64 = branch["g_fr"]
        b_fr::Float64 = branch["b_fr"]
        g_to::Float64 = branch["g_to"]
        b_to::Float64 = branch["b_to"]

        # The formula below match the PM implementation,
        #  given Ohm's in the form:
        #    gff * wf + gft * wr + bft * wi == pf
        #   -bff * wf - bft * wr + gft * wi == qf
        #    gtt * wt + gtf * wr - btf * wi == pt
        #   -btt * wt - btf * wr - gtf * wi == qt
        # Note: OPOMO and Gurobi optimods use the same form of Ohm's law
        gff[e] = (g+g_fr) / ttm
        gft[e] = (-g*tr+b*ti) / ttm
        gtf[e] = (-g*tr-b*ti) / ttm
        gtt[e] = (g+g_to)
        bff[e] = (b+b_fr) / ttm
        bft[e] = (-b*tr-g*ti) / ttm
        btf[e] = (-b*tr+g*ti) / ttm
        btt[e] = (b+b_to)

        # Angle deviation constraints
        dvamin[e] = branch["angmin"]
        dvamax[e] = branch["angmax"]

        # Thermal limits
        smax[e] = branch["rate_a"]

        branch_status[e] = branch["br_status"] == 1

        # Branch incidence matrix
        A_i[e] = e
        A_j[e] = i
        A_v[e] = 1

        A_i[E+e] = e
        A_j[E+e] = j
        A_v[E+e] = -1
    end
    sort!.(bus_arcs_fr)
    sort!.(bus_arcs_to)

    A = sparse(A_i, A_j, A_v, E, N)

    return OPFData(
        network["name"], network["baseMVA"],
        N, E, G, L, A, Ag,
        vnom, vmin, vmax, gs, bs, pd, qd,
        bus_arcs_fr, bus_arcs_to, bus_gens, ref_bus,
        pgmin, pgmax,
        qgmin, qgmax,
        c0, c1, c2,
        gen_status,
        bus_fr, bus_to,
        branch_g, branch_b,
        gff, gft, gtf, gtt,
        bff, bft, btf, btt,
        smax, dvamin, dvamax,
        branch_status,
    )
end

function to_dict(data::OPFData)
    d = Dict{String,Any}()
    for field in fieldnames(OPFData)
        v = getfield(data, field)
        if isa(v, SparseMatrixCSC)
            I, J, V = findnz(v)
            M, N = size(v)
            d[string(field)] = Dict("I" => I, "J" => J, "V" => V, "M" => M, "N" => N)
        elseif isa(v, AbstractArray)
            d[string(field)] = copy(v)
        else
            d[string(field)] = v
        end
    end
    return d
end

include("utils.jl")
include("acp.jl")      # ACPPowerModel
include("dcp.jl")      # DCPPowerModel
include("ed.jl")       # EconomicDispatch
include("socwr.jl")    # SOCWRPowerModel & SOCWRConicPowerModel

# Contains a list of all supported OPF models
const SUPPORTED_OPF_MODELS = Type{<:Formulation}[
    ACOPF,
    DCOPF,
    EconomicDispatch,
    SOCOPFQuad,
    PowerModels.SOCWRConicPowerModel,
]

# A name --> type dictionary
# Used for passing the OPF type as a string (e.g. through config file)
const OPF2TYPE = Dict{String,Type{<:Formulation}}(
    "ACOPF" => ACOPF,
    "DCOPF" => DCOPF,
    "EconomicDispatch" => EconomicDispatch,
    "SOCWRPowerModel" => SOCOPFQuad,
    "SOCWRConicPowerModel" => PowerModels.SOCWRConicPowerModel,
)

function solve!(opf::OPFModel{<:Formulation})
    optimize!(opf.model; _differentiation_backend = MathOptSymbolicAD.DefaultBackend())
end
