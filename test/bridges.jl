function test_bridges_path(n)
    data = Dict{String,Any}(
        "basic_network" => true,
    )
    # Create `n` buses
    data["bus"] = Dict{String,Any}(
        "$i" => Dict{String,Any}()
        for i in 1:n
    )

    # Create `n-1` branches (i, i+1)
    # The graph is a path [1]--[2]--...--[n] and every branch is a bridge
    data["branch"] = Dict{String,Any}(
        "$i" => Dict{String,Any}(
            "f_bus" => i,
            "t_bus" => i+1,
        )
        for i in 1:(n-1)
    )
    b = OPFGenerator.bridges(data)
    @test all(values(b))

    # Add one branch between `n` and `1`
    # --> the graph is a cycle and no branch is a bridge
    data["branch"]["$n"] = Dict{String,Any}(
        "f_bus" => 1,
        "t_bus" => n,
    )
    b = OPFGenerator.bridges(data)
    @test all(!b["$i"] for i in 1:n)

    return nothing
end

function test_bridge_star(n)
    data = Dict{String,Any}(
        "basic_network" => true,
    )
    # Create `n` buses
    data["bus"] = Dict{String,Any}(
        "$i" => Dict{String,Any}()
        for i in 1:n+1
    )
    # Create `n` branches between bus `n+1` and bus `i`
    data["branch"] = Dict(
        "$i" => Dict(
            "f_bus" => i,
            "t_bus" => n+1,
        )
        for i in 1:n
    )
    b = OPFGenerator.bridges(data)
    @test all(values(b))

    # Add branch between bus `i` and bus `n+1`
    # --> the corresponding branch should no longer be a bridge
    data["branch"]["$(n+1)"] = Dict(
        "f_bus" => n+1,
        "t_bus" => n+1
    )
    for i in 1:n
        data["branch"]["$(n+1)"]["t_bus"] = i
        b = OPFGenerator.bridges(data)
        @test !b["$i"]
        @test !b["$(n+1)"]
        @test all(b["$j"] for j in 1:n if j != i)
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
    data = Dict{String,Any}("basic_network" => false)
    @test_throws ErrorException OPFGenerator.bridges(data)

    data["basic_network"] => true
    data["bus"] = Dict("1" => Dict(), "2" => Dict())
    data["branch"] = Dict("1" => Dict("f_bus" => 1, "t_bus" => 1))
    @test_throws ErrorException OPFGenerator.bridges(data)

    return nothing
end

@testset "Bridges" begin
    @testset test_bridges_path(4)
    @testset test_bridge_star(3)
    @testset test_bridge_edgecases()
end
