using OPFGenerator: LazyPTDF, FullPTDF, OPFData, compute_flow!, ptdf_row

function test_ptdf(network, data)
    N = length(network["bus"])
    E = length(network["branch"])
    p = ones(N)
    f = zeros(E)

    Φ_pm = PowerModels.calc_basic_ptdf_matrix(network)
    fpm = Φ_pm * p

    for PTDF in [FullPTDF, LazyPTDF]
        Φ = PTDF(data)
        f = zeros(E)
        compute_flow!(f, p, Φ)
        @test isapprox(f, fpm; atol=1e-6)

        # Test ptdf_row
        for e in 1:E
            row = ptdf_row(Φ, e)
            @test row == Φ.matrix[e, :]
            @test isapprox(row, Φ_pm[e, :]; atol=1e-6)
        end
    end

    return nothing
end
