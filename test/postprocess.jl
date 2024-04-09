function test_tensorize()
    # Vector of scalars
    V = [1.0, 2.0, 3.0, 4.0]
    @test OPFGenerator.tensorize(V) == [1.0 2.0 3.0 4.0]
    V = ["1.0", "2.0", "3.0", "4.0", "e", "f"]
    @test OPFGenerator.tensorize(V) == ["1.0" "2.0" "3.0" "4.0" "e" "f"]
    
    # Vector of vectors
    V = [(10*i) .+ collect(1:4) for i in 1:3]
    @test OPFGenerator.tensorize(V) == [
        10 21 31;
        11 22 32;
        13 23 33;
        14 24 34;
    ]

    # Vector of matrices
    V = [(i .* ones(2, 3)) for i in 1:4]
    M = OPFGenerator.tensorize(V)
    @test size(M) == (2, 3, 4)
    @test all(M[:, :, i] == V[i] for i in 1:4)

    return nothing
end

@testset test_tensorize()
