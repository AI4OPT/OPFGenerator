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
