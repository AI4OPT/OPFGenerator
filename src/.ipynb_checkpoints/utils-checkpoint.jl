function _set_loads!(data, pd, qd)
    L = length(data["load"])
    length(pd) == L || throw(DimensionMismatch())
    length(qd) == L || throw(DimensionMismatch())
    
    for i in 1:L
        ldat = data["load"]["$i"]
        ldat["pd"] = pd[i]
        ldat["qd"] = qd[i]
    end

    return nothing
end