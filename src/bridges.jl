using Graphs

"""
    bridges(data)

Identify whether each branch is a bridge.

The input data must be in basic format.

A branch is a bridge if removing it renders the network disconnected.
Returns a dictionary `res::Dict{String,Bool}` such that
    `res[br]` is true if branch `br` is a bridge, and `false` otherwise.
"""
function bridges(data::Dict)
    get(data, "basic_network", false) || error(
        """Invalid data: network data must be in basic format.
        Call `make_basic_network(data)` before calling this function"""
    )

    N::Int = length(data["bus"])
    E::Int = length(data["branch"])

    # Build the adjacency graph...
    G = Graph(N)
    # ... and record which branches have dupplicates
    D = Dict{Tuple{Int,Int},Vector{Int}}()
    for e in 1:E
        br::Dict = data["branch"]["$e"]
        i::Int = br["f_bus"]
        j::Int = br["t_bus"]

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
    return Dict("$e" => res[e] for e in 1:E)
end
