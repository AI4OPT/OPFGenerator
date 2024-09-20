function voltage_phasor_bounds(data::OPFData)
    E = data.E
    vmin, vmax = data.vmin, data.vmax
    dvamin, dvamax = data.dvamin, data.dvamax
    bus_fr, bus_to = data.bus_fr, data.bus_to
    return voltage_phasor_bounds(E, vmin, vmax, dvamin, dvamax, bus_fr, bus_to)
end

function voltage_phasor_bounds(E, vmin, vmax, dvamin, dvamax, bus_fr, bus_to)
    wr_min = zeros(Float64, E)
    wr_max = zeros(Float64, E)
    wi_min = zeros(Float64, E)
    wi_max = zeros(Float64, E)
    for e in 1:E
        i = bus_fr[e]
        j = bus_to[e]
        cosmin = cos(dvamin)
        cosmax = cos(dvamax)
        sinmin = sin(dvamin)
        sinmax = sin(dvamax)

        wr_min[e] = (
            dvamin >= 0 ? vmin[i] * vmin[j] * cosmax :
            dvamax <= 0 ? vmin[i] * vmin[j] * cosmin :
            vmin[i] * vmin[j] * min(cosmin, cosmax)
        )
        wr_max[e] = (
            dvamin >= 0 ? vmax[i] * vmax[j] * cosmin :
            dvamax <= 0 ? vmax[i] * vmax[j] * cosmax :
            vmax[i] * vmax[j] * 1.0
        )
        wi_min[e] = (
            dvamin >= 0 ? vmin[i] * vmin[j] * sinmin :
            dvamax <= 0 ? vmax[i] * vmax[j] * sinmax :
            vmax[i] * vmax[j] * sinmin
        )
        wi_max[e] = (
            dvamin >= 0 ? vmax[i] * vmax[j] * sinmax :
            dvamax <= 0 ? vmin[i] * vmin[j] * sinmax :
            vmax[i] * vmax[j] * sinmax
        )
    end

    return wr_min, wr_max, wi_min, wi_max
end
