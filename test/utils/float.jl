function test_float_conversion()
    D = Dict(
        "a" => 1,
        "b" => 2.0,
        "c" => [3.1, 3.2, 3.3],
        "d" => Dict(
            "d1" => ones(BigFloat, (1,)),
            "d2" => ones(BigFloat, (1,2)),
            "d3" => ones(BigFloat, (3,1,2)),
        ),
        "e" => "hello, world!",
        "f" => [1, 2, 3],
    )

    # Test Float32 conversion
    # Use `===` to check value _and_ type when applicable,
    #   and to check that we didn't copy any data that wasn't converted
    d32 = OPFGenerator.convert_float_data(D, Float32)
    @test isa(d32, Dict{String,Any})
    @test d32["a"] === 1
    @test d["b"] === 2f0
    @test isa(d32["c"], Array{Float32,1})
    @test d32["c"] == [3.1f0, 3.2f0, 3.3f0]
    @test isa(d32["d"], Dict{String,Any})
    @test length(d32["d"]) == 3
    @test isa(d32["d"]["d1"], Array{Float32,1})
    @test d32["d"]["d1"] == ones(Float32, (1,))
    @test isa(d32["d"]["d2"], Array{Float32,2})
    @test d32["d"]["d1"] == ones(Float32, (1,2))
    @test isa(d32["d"]["d3"], Array{Float32,3})
    @test d32["d"]["d1"] == ones(Float32, (3,1,2))
    @test d32["e"] === D["e"]
    @test d32["f"] === D["f"]

    # Test Float64 conversion
    d64 = OPFGenerator.convert_float_data(D, Float64)
    @test isa(d64, Dict{String,Any})
    @test d64["a"] === 1
    @test d64["b"] === 2.0
    @test d64["c"] === D["c"]
    @test isa(d64["d"], Dict{String,Any})
    @test length(d64["d"]) == 3
    @test isa(d64["d"]["d1"], Array{Float64,1})
    @test d64["d"]["d1"] == ones(Float64, (1,))
    @test isa(d64["d"]["d2"], Array{Float64,2})
    @test d64["d"]["d1"] == ones(Float64, (1,2))
    @test isa(d64["d"]["d3"], Array{Float64,3})
    @test d64["d"]["d1"] == ones(Float64, (3,1,2))
    @test d64["e"] === D["e"]
    @test d64["f"] === D["f"]

    # Test argument sanity checks
    @test_throws ErrorException OPFGenerator.convert_float_data(Dict{Int,Any}(), Float32)
    @test_throws ErrorException OPFGenerator.convert_float_data(Dict{String,Any}(), Int)

    return nothing
end

@testset "Float conversion" begin
    test_float_conversion()
end
