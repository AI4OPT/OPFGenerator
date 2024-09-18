## arcs

ArcType = Tuple{Int, Int, Int}

function make_arcs_from(data::Dict{String,Any})
    return [
        (l, data["branch"]["$l"]["f_bus"], data["branch"]["$l"]["t_bus"])
        for l in 1:length(data["branch"])
    ]
end

function make_arcs_to(data::Dict{String,Any})
    return [
        (l, data["branch"]["$l"]["t_bus"], data["branch"]["$l"]["f_bus"])
        for l in 1:length(data["branch"])
    ]
end

function make_arcs_to(arcs_from::Vector{ArcType})
    return [(l, j, i) for (l, i, j) in arcs_from]
end

function make_arcs(data::Dict{String,Any})
    arcs_from = make_arcs_from(data)
    return make_arcs(arcs_from)
end

function make_arcs(arcs_from::Vector{ArcType})
    return [arcs_from; make_arcs_to(arcs_from)]
end

## bus_*

function make_bus_loads(data::Dict{String,Any})
    bus_loads = Vector{Vector{Int}}(undef, length(data["bus"]))
    L = length(data["load"])
    for b in 1:length(data["bus"])
        bus_loads[b] = [
            l for l in 1:L if data["load"]["$l"]["load_bus"] == b
        ]
    end
    return bus_loads
end

function make_bus_shunts(data::Dict{String,Any}) # TODO: vector of vectors
    bus_shunts = Vector{Vector{Int}}(undef, length(data["bus"]))
    S = length(data["shunt"])
    for b in 1:length(data["bus"])
        bus_shunts[b] = [
            s for s in 1:S if data["shunt"]["$s"]["shunt_bus"] == b
        ]
    end
    return bus_shunts
end

function make_bus_gens(data::Dict{String,Any})
    bus_gens = Vector{Vector{Int}}(undef, length(data["bus"]))
    G = length(data["gen"])
    for b in 1:length(data["bus"])
        bus_gens[b] = [
            g for g in 1:G if data["gen"]["$g"]["gen_bus"] == b
        ]
    end
    return bus_gens
end

function make_bus_arcs(data::Dict{String,Any})
    arcs_from = make_arcs_from(data)
    return make_bus_arcs(data, arcs_from)
end

function make_bus_arcs(data::Dict{String,Any}, arcs_from::Vector{ArcType})
    bus_arcs = Vector{Vector{ArcType}}(undef, length(data["bus"]))
    for b in 1:length(data["bus"])
        bus_arcs[b] = [
            # if from bus, add arc as is
            # if to bus, swap from and to
            i == b ? (l, i, j) : (l, j, i)
            for (l, i, j) in arcs_from
            if (i == b || j == b)
        ]
    end
    return bus_arcs
end

function make_ref_buses(data::Dict{String,Any})
    return [
        b for b in 1:length(data["bus"])
        if data["bus"]["$b"]["bus_type"] == 3
    ]
end

## active_*

function make_active_branches(data::Dict{String,Any})
    return [
        l for l in 1:length(data["branch"])
        if data["branch"]["$l"]["br_status"] == 1
    ]
end

function make_active_generators(data::Dict{String,Any})
    return [
        g for g in 1:length(data["gen"])
        if data["gen"]["$g"]["gen_status"] == 1
    ]
end

""" 
    fix_variables(status, vars)

Fix variables by setting their lower and upper bounds to zero if the corresponding status is zero.
Note this differs from JuMP.fix which creates an equality constraint and removes bounds.
"""
function zero_variables(status::Vector{Int}, vars)
    @assert length(status) == length(vars) "status and vars must have the same length"
    for (s, v) in zip(status, vars)
        if s == 0
            JuMP.set_lower_bound(v, 0.0)
            JuMP.set_upper_bound(v, 0.0)
        end
    end
end


function calc_buspair_voltage_bounds!(data::Dict{String,Any})
    E = length(data["branch"])
    
    data["buspair"] = buspairs = Dict{Tuple{Int,Int},Dict{String,Float64}}()
    for e in 1:E
        i = data["branch"]["$e"]["f_bus"]
        j = data["branch"]["$e"]["t_bus"]

        if !haskey(buspairs, (i,j))
            buspairs[(i,j)] = Dict{String,Float64}()
            buspairs[(i,j)]["angmin"] = data["branch"]["$e"]["angmin"]
            buspairs[(i,j)]["angmax"] = data["branch"]["$e"]["angmax"]
        else
            buspairs[(i,j)]["angmin"] = max(buspairs[(i,j)]["angmin"], data["branch"]["$e"]["angmin"])
            buspairs[(i,j)]["angmax"] = min(buspairs[(i,j)]["angmax"], data["branch"]["$e"]["angmax"])
        end

    end

    for (i,j) in keys(buspairs)
        buspair = buspairs[(i,j)]

        busi = data["bus"]["$i"]
        busj = data["bus"]["$j"]

        if buspair["angmin"] >= 0
            buspairs[(i,j)]["wr_min"] = busi["vmin"]*busj["vmin"]*cos(buspair["angmax"])
            buspairs[(i,j)]["wr_max"] = busi["vmax"]*busj["vmax"]*cos(buspair["angmin"])
            buspairs[(i,j)]["wi_min"] = busi["vmin"]*busj["vmin"]*sin(buspair["angmin"])
            buspairs[(i,j)]["wi_max"] = busi["vmax"]*busj["vmax"]*sin(buspair["angmax"])
        end
        if buspair["angmax"] <= 0
            buspairs[(i,j)]["wr_min"] = busi["vmin"]*busj["vmin"]*cos(buspair["angmin"])
            buspairs[(i,j)]["wr_max"] = busi["vmax"]*busj["vmax"]*cos(buspair["angmax"])
            buspairs[(i,j)]["wi_min"] = busi["vmax"]*busj["vmax"]*sin(buspair["angmin"])
            buspairs[(i,j)]["wi_max"] = busi["vmin"]*busj["vmin"]*sin(buspair["angmax"])
        end
        if buspair["angmin"] < 0 && buspair["angmax"] > 0
            buspairs[(i,j)]["wr_min"] = busi["vmin"]*busj["vmin"]*min(cos(buspair["angmin"]), cos(buspair["angmax"]))
            buspairs[(i,j)]["wr_max"] = busi["vmax"]*busj["vmax"]*1.0
            buspairs[(i,j)]["wi_min"] = busi["vmax"]*busj["vmax"]*sin(buspair["angmin"])
            buspairs[(i,j)]["wi_max"] = busi["vmax"]*busj["vmax"]*sin(buspair["angmax"])
        end
    end
end