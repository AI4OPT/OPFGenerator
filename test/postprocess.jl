function test_tensorize()
    # Vector of scalars
    V = [1.0, 2.0, 3.0, 4.0]
    @test OPFGenerator.tensorize(V) == [1.0, 2.0, 3.0, 4.0]
    V = ["1.0", "2.0", "3.0", "4.0", "e", "f"]
    @test OPFGenerator.tensorize(V) == ["1.0", "2.0", "3.0", "4.0", "e", "f"]
    
    # Vector of vectors
    V = [(10*i) .+ collect(1:4) for i in 1:3]
    @test OPFGenerator.tensorize(V) == [
        11 21 31;
        12 22 32;
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

function test_merge_h5_array()
    # Cannot merge empty collection
    V = Array{Float64,1}[]
    @test_throws ErrorException OPFGenerator._merge_h5(V)

    # Check that trying to merge arrays of wrong sizes yields an error
    V = [ones(1, 2, 3), ones(2, 1, 3)]
    @test_throws DimensionMismatch OPFGenerator._merge_h5(V)

    # Merge arrays of multiple dimensions
    for N in 1:4
        # ⚠ the memory below is `N!`, so keep `N` under control (no more than 6)
        V = [i .* ones(Float64, 1:N...) for i in 1:4]

        # Check for type stability
        M = @inferred Array{Float64,N} OPFGenerator._merge_h5(V)

        # Checck output
        @test size(M) == (1:(N-1)..., 4*N)
        # `_merge_h5` should concatenate the arrays along the last dimension
        # (but uses a more efficient implementation for dealing with many arrays)
        @test M == cat(V...; dims=ndims(M))
    end

    # merge with different minibatch sizes
    # a `Base.stack`-based implementation should fail this
    V = [
        [100 * i1 + 10*i2 + i3 for i1 in 1:2, i2 in 1:3, i3 in 1:2],
        [100 * i1 + 10*i2 + (2+i3) for i1 in 1:2, i2 in 1:3, i3 in 1:4],
        [100 * i1 + 10*i2 + (6+i3) for i1 in 1:2, i2 in 1:3, i3 in 1:3],
    ]
    M = OPFGenerator._merge_h5(V)
    @test M == [100 * i1 + 10*i2 + i3 for i1 in 1:2, i2 in 1:3, i3 in 1:9]

    return nothing
end

@testset test_tensorize()
@testset test_merge_h5_array()
