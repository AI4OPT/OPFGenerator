struct OPFData
    N::Int  # number of buses
    E::Int  # number of branches
    G::Int  # number of generators

    # Bus data
    vmin::Vector{Float64}
    vmax::Vector{Float64}
    gs::Vector{Float64}  # shunt
    bs::Vector{Float64}  # shunt
    pd::Vector{Float64}  # nodal active power load
    qd::Vector{Float64}  # nodal reactive power load
    bus_arcs_fr::Vector{Vector{Int}}  # indices of branches exiting the bus
    bus_arcs_to::Vector{Vector{Int}}  # indices of branches entering the bus
    bus_gens::Vector{Vector{Int}}  # indices of generators at the bus
    ref_bus::Int  # reference bus

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
    @assert network["per_unit"] == true "Network data is not per-unit scaled."
    @assert network["baseMVA"] == 100.0 "Base MVA is not 100.0."

    N = length(network["bus"])
    E = length(network["branch"])
    G = length(network["gen"])

    # Bus data
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
    end
    # sort everything again
    sort!.(bus_gens)

    # Branch data
    bus_fr = zeros(Int, E)
    bus_to = zeros(Int, E)
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
    end
    sort!.(bus_arcs_fr)
    sort!.(bus_arcs_to)

    return OPFData(
        N, E, G,
        vmin, vmax, gs, bs, pd, qd,
        bus_arcs_fr, bus_arcs_to, bus_gens, ref_bus,
        pgmin, pgmax,
        qgmin, qgmax,
        c0, c1, c2,
        gen_status,
        bus_fr, bus_to,
        gff, gft, gtf, gtt,
        bff, bft, btf, btt,
        smax, dvamin, dvamax,
        branch_status,
    )
end

function buspair_voltage_bounds(data::OPFData)
    E = data.E
    vmin, vmax = data.vmin, data.vmax
    dvamin, dvamax = data.dvamin, data.dvamax
    bus_fr, bus_to = data.bus_fr, data.bus_to
    
    buspairs = Dict{Tuple{Int,Int},Dict{String,Any}}()
    for e in 1:E
        i = bus_fr[e]
        j = bus_to[e]

        if !haskey(buspairs, (i,j))
            buspairs[(i,j)] = buspair = Dict{String,Any}()
            buspair["angmin"] = dvamin[e]
            buspair["angmax"] = dvamax[e]
            buspair["edges"] = [e]
        else
            buspair = buspairs[(i,j)]
            buspair["angmin"] = max(buspair["angmin"], dvamin[e])
            buspair["angmax"] = min(buspair["angmax"], dvamax[e])
            push!(buspair["edges"], e)
        end
    end

    wr_min = zeros(Float64, E)
    wr_max = zeros(Float64, E)
    wi_min = zeros(Float64, E)
    wi_max = zeros(Float64, E)

    for (i,j) in keys(buspairs)
        buspair = buspairs[(i,j)]
        cosmin = cos(buspair["angmin"])
        cosmax = cos(buspair["angmax"])
        sinmin = sin(buspair["angmin"])
        sinmax = sin(buspair["angmax"])

        bp_wr_min = (
            buspair["angmin"] >= 0 ? vmin[i] * vmin[j] * cosmax :
            buspair["angmax"] <= 0 ? vmin[i] * vmin[j] * cosmin :
            vmin[i] * vmin[j] * min(cosmin, cosmax)
        )
        bp_wr_max = (
            buspair["angmin"] >= 0 ? vmax[i] * vmax[j] * cosmin :
            buspair["angmax"] <= 0 ? vmax[i] * vmax[j] * cosmax :
            vmax[i] * vmax[j] * 1.0
        )
        bp_wi_min = (
            buspair["angmin"] >= 0 ? vmin[i] * vmin[j] * sinmin :
            buspair["angmax"] <= 0 ? vmax[i] * vmax[j] * sinmax :
            vmax[i] * vmax[j] * sinmin
        )
        bp_wi_max = (
            buspair["angmin"] >= 0 ? vmax[i] * vmax[j] * sinmax :
            buspair["angmax"] <= 0 ? vmin[i] * vmin[j] * sinmax :
            vmax[i] * vmax[j] * sinmax
        )

        for e in buspair["edges"]
            wr_min[e] = bp_wr_min
            wr_max[e] = bp_wr_max
            wi_min[e] = bp_wi_min
            wi_max[e] = bp_wi_max
        end
    end

    return wr_min, wr_max, wi_min, wi_max
end
