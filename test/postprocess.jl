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

function test_dedupe_sort_h5()
    # Sample data for testing
    D = Dict{String,Any}(
        "meta" => Dict("seed" => [1, 1, 4, 3, 3, 2]),
        # ⚠ the values below are intentionally slightly different across
        #   repetition of the same random seed; this is for testing purposes
        "primal" => Dict(
            "pg" => [10.1, 10.2, 40.1, 30.1, 30.2, 20.1],
            "qg" => [11.1, 11.2, 41.1, 31.1, 31.2, 21.1],
        ),
        "dual" => Dict(
            "lambda" => [
                101.1 101.2 401.1 301.1 301.2 201.1;
                102.1 102.2 402.1 302.1 302.2 202.1;
            ],
            "mu" => [
                111.1 111.2 411.1 311.1 311.2 211.1;
                112.1 112.2 412.1 312.1 312.2 212.1;
            ],
        )
    )

    # dedupe
    @testset "dedupe" begin
    D_ = OPFGenerator._dedupe_h5!(deepcopy(D))
        @test D_["meta"]["seed"] == [1, 4, 3, 2]
        @test D_["primal"]["pg"] == [10.1, 40.1, 30.1, 20.1]
        @test D_["primal"]["qg"] == [11.1, 41.1, 31.1, 21.1]
        @test D_["dual"]["lambda"] == [
            101.1 401.1 301.1 201.1;
            102.1 402.1 302.1 202.1;
        ]
        @test D_["dual"]["mu"] == [
            111.1 411.1 311.1 211.1;
            112.1 412.1 312.1 212.1;
        ]
    end

    # sort
    @testset "sort" begin
        D_ = OPFGenerator._sort_h5!(deepcopy(D))
        @test D_["meta"]["seed"] == [1, 1, 2, 3, 3, 4]
        # we cannot guarantee that the order of repeated seeds is preserved,
        #   so we only check the integer parts 
        #   (decimals only encode the repetition of the seeds)
        @test floor.(D_["primal"]["pg"]) == [10, 10, 20, 30, 30, 40]
        @test floor.(D_["primal"]["qg"]) == [11, 11, 21, 31, 31, 41]
        @test floor.(D_["dual"]["lambda"]) == [
            101 101 201 301 301 401;
            102 102 202 302 302 402;
        ]
        @test floor.(D_["dual"]["mu"]) == [
            111 111 211 311 311 411;
            112 112 212 312 312 412;
        ]
    end
    
    # dedupe and sort
    @testset "dedupe and sort" begin
        D_ = OPFGenerator._dedupe_and_sort_h5!(deepcopy(D))
        # the dedupe part should keep only the first occurence of each seed,
        #   so here we do test with decimals values
        @test D_["meta"]["seed"] == [1, 2, 3, 4]
        @test D_["primal"]["pg"] == [10.1, 20.1, 30.1, 40.1]
        @test D_["primal"]["qg"] == [11.1, 21.1, 31.1, 41.1]
        @test D_["dual"]["lambda"] == [
            101.1 201.1 301.1 401.1;
            102.1 202.1 302.1 402.1;
        ]
        @test D_["dual"]["mu"] == [
            111.1 211.1 311.1 411.1;
            112.1 212.1 312.1 412.1;
        ]
    end

    return nothing
end

@testset test_tensorize()
@testset test_merge_h5_array()
@testset test_dedupe_sort_h5()
