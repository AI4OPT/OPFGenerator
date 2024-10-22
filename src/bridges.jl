using Graphs

"""
    bridges(data)

Identify whether each branch is a bridge.

The input data must be in basic format.

A branch is a bridge if removing it renders the network disconnected.
Returns a dictionary `res::Dict{String,Bool}` such that
    `res[br]` is true if branch `br` is a bridge, and `false` otherwise.
"""
function bridges(data::OPFData)
    N, E = data.N, data.E
    bus_fr, bus_to = data.bus_fr, data.bus_to

    return bridges(N, E, bus_fr, bus_to)
end

function bridges(N::Int, E::Int, bus_fr::Vector{Int}, bus_to::Vector{Int})
    # Build the adjacency graph...
    G = Graph(N)
    # ... and record which branches have dupplicates
    D = Dict{Tuple{Int,Int},Vector{Int}}()
    for e in 1:E
        i = bus_fr[e]
        j = bus_to[e]

        # Ensure all edges are of the form i --> j with i ≤ j
        _i, _j = extrema((i, j))
        (1 <= _i < _j <= N) || error("Invalid branch #$e=($i, $j)")
        # Check if there's a dupplicate
        u = get!(D, (_i, _j), Int[])
        push!(u, e)
        add_edge!(G, _i, _j)
    end

    # Compute bridges
    B = Graphs.bridges(G)

    # A branch is a bridge if:
    # * It does not have a dupplicate edge, and
    # * It is a bridge in the subgraph with no multi-edges
    res = falses(E)
    for edge in B
        i, j = edge.src, edge.dst
        u = D[i, j]
        if length(u) == 1
            res[u] .= true
        end
    end
    return res
end