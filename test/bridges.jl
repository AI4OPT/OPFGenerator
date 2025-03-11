function test_bridges_path(n)
    # Create `n-1` branches (i, i+1)
    # The graph is a path [1]--[2]--...--[n] and every branch is a bridge
    bus_fr = [i for i in 1:n-1]
    bus_to = [i+1 for i in 1:n-1]
    b = PGLearn.bridges(n, n-1, bus_fr, bus_to)
    @test all(values(b))

    # Add one branch between `n` and `1`
    # --> the graph is a cycle and no branch is a bridge
    append!(bus_fr, [n])
    append!(bus_to, [1])
    b = PGLearn.bridges(n, n, bus_fr, bus_to)
    @test all(!b[i] for i in 1:n)

    return nothing
end

function test_bridge_star(n)
    bus_fr = [i for i in 1:n]
    bus_to = [n+1 for _ in 1:n]
    b = PGLearn.bridges(n+1, n, bus_fr, bus_to)
    @test all(values(b))

    # Add branch between bus `i` and bus `n+1`
    # --> the corresponding branch should no longer be a bridge
    append!(bus_fr, [n+1])
    append!(bus_to, [n+1])
    for i in 1:n
        bus_to[n+1] = i
        b = PGLearn.bridges(n+1, n+1, bus_fr, bus_to)
        @test !b[i]
        @test !b[n+1]
        @test all(b[j] for j in 1:n if j != i)
    end

    return nothing
end

"""
    test_bridge_edgecases()

Test following edge-cases in the `bridges` function:
* data is not in basic format
* loops, i.e., branch of the form (i, i)
"""
function test_bridge_edgecases()
    @test_throws ErrorException PGLearn.bridges(2, 1, [1], [1])

    return nothing
end

@testset "Bridges" begin
    @testset test_bridges_path(4)
    @testset test_bridge_star(3)
    @testset test_bridge_edgecases()
end