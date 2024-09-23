using OPFGenerator: LazyPTDF, FullPTDF, OPFData, compute_flow!, ptdf_row

function test_ptdf_full(network, data)
    N = length(network["bus"])
    E = length(network["branch"])
    Φ_pm = PowerModels.calc_basic_ptdf_matrix(network)
    p = randn(N)
    f = zeros(E)

    Φ = FullPTDF(data)
    fpm = Φ_pm * p
    compute_flow!(f, p, Φ)
    @test isapprox(f, fpm; atol=1e-6)

    # Test ptdf_row
    for e in 1:E
        row = ptdf_row(Φ, e)
        @test row == Φ.matrix[e, :]
        @test isapprox(row, Φ_pm[e, :]; atol=1e-6)
    end

    return nothing
end

function test_ptdf_lazy(network, data)
    N = length(network["bus"])
    E = length(network["branch"])

    p = randn(N)
    f = zeros(E)
    Φ_pm = PowerModels.calc_basic_ptdf_matrix(network)
    fpm = Φ_pm * p
    
    @testset "$solver" for solver in [:ldlt, :cholesky]
        Φ = LazyPTDF(data; solver=solver)

        # Check power flow computation
        compute_flow!(f, p, Φ)
        @test isapprox(f, fpm; atol=1e-6)

        # Test ptdf_row
        for e in 1:E
            row = ptdf_row(Φ, e)
            @test isapprox(row, Φ_pm[e, :]; atol=1e-6)
        end
    end

    return nothing
end