function voltage_phasor_bounds(data::OPFData)
    E = data.E
    vmin, vmax = data.vmin, data.vmax
    dvamin, dvamax = data.dvamin, data.dvamax
    bus_fr, bus_to = data.bus_fr, data.bus_to
    vfmin = [vmin[bus_fr[e]] for e in 1:E]
    vfmax = [vmax[bus_fr[e]] for e in 1:E]
    vtmin = [vmin[bus_to[e]] for e in 1:E]
    vtmax = [vmax[bus_to[e]] for e in 1:E]
    return voltage_phasor_bounds(E, vfmin, vfmax, vtmin, vtmax, dvamin, dvamax)
end

function voltage_phasor_bounds(E, vfmin, vfmax, vtmin, vtmax, dvamin, dvamax)
    wr_min = zeros(Float64, E)
    wr_max = zeros(Float64, E)
    wi_min = zeros(Float64, E)
    wi_max = zeros(Float64, E)

    cosmin = cos.(dvamin)
    cosmax = cos.(dvamax)
    sinmin = sin.(dvamin)
    sinmax = sin.(dvamax)

    for e in 1:E
        wr_min[e] = (
            dvamin[e] >= 0 ? vfmin[e] * vtmin[e] * cosmax[e] :
            dvamax[e] <= 0 ? vfmin[e] * vtmin[e] * cosmin[e] :
            vfmin[e] * vtmin[e] * min(cosmin[e], cosmax[e])
        )
        wr_max[e] = (
            dvamin[e] >= 0 ? vfmax[e] * vtmax[e] * cosmin[e] :
            dvamax[e] <= 0 ? vfmax[e] * vtmax[e] * cosmax[e] :
            vfmax[e] * vtmax[e] * 1.0
        )
        wi_min[e] = (
            dvamin[e] >= 0 ? vfmin[e] * vtmin[e] * sinmin[e] :
            dvamax[e] <= 0 ? vfmax[e] * vtmax[e] * sinmax[e] :
            vfmax[e] * vtmax[e] * sinmin[e]
        )
        wi_max[e] = (
            dvamin[e] >= 0 ? vfmax[e] * vtmax[e] * sinmax[e] :
            dvamax[e] <= 0 ? vfmin[e] * vtmin[e] * sinmax[e] :
            vfmax[e] * vtmax[e] * sinmax[e]
        )
    end

    return wr_min, wr_max, wi_min, wi_max
end
