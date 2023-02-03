using Graphs

function bridges(data::Dict)
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
