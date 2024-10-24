using Distributions
using LinearAlgebra
using TOML

function test_glocal()
    d = OPFGenerator.Glocal(
        Uniform(0.0, 1.0),
        Distributions.MvNormal(zeros(4), Diagonal(ones(4)))
    )

    @test eltype(d) == Float64
    @test length(d) == 4

    @test size(rand(d)) == (4,)
    @test size(rand(d, 1)) == (4, 1)
    @test size(rand(d, 2)) == (4, 2)

    return nothing
end

function test_ScaledLogNormal()
    d = ScaledLogNormal(0.8, 1.2, 0.05 .* ones(3))

    @test length(d) == 3

    @test isa(d, OPFGenerator.Glocal)
    @test d.d_α == Uniform(0.8, 1.2)
    @test isa(d.d_η, Distributions.MvLogNormal)

    # Sanity checks
    @test_throws ErrorException ScaledLogNormal(0.8, 0.7, ones(3))   # l > u
    @test_throws ErrorException ScaledLogNormal(0.8, 1.2, -ones(3))  # σ < 0

    return nothing
end

function test_ScaledUniform()
    d = ScaledUniform(0.8, 1.2, 0.05 .* ones(5))

    @test length(d) == 5

    @test isa(d, OPFGenerator.Glocal)
    @test d.d_α == Uniform(0.8, 1.2)
    @test isa(d.d_η, Distributions.Product)

    # Sanity checks
    @test_throws ErrorException ScaledUniform(0.8, 0.7, ones(3))   # l > u
    @test_throws ErrorException ScaledUniform(0.8, 1.2, -ones(3))  # σ < 0

    return nothing
end

function test_LoadScaler()
    data = OPFGenerator.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))

    # ScaledLogNormal
    options = Dict(
        "noise_type" => "ScaledLogNormal",
        "l" => 0.8,
        "u" => 1.2,
        "sigma" => 0.05,
    )
    ls = LoadScaler(data, options)
    @test ls.d.d_α == Uniform(0.8, 1.2)
    @test isa(ls.d.d_η, MvLogNormal)
    @test ls.pd_ref == data.pd
    @test ls.qd_ref == data.qd

    # ScaledUniform
    options = Dict(
        "noise_type" => "ScaledUniform",
        "l" => 0.7,
        "u" => 1.5,
        "sigma" => 0.05,
    )
    ls = LoadScaler(data, options)
    @test ls.d.d_α == Uniform(0.7, 1.5)
    @test isa(ls.d.d_η, Distributions.Product)
    @test ls.pd_ref == data.pd
    @test ls.qd_ref == data.qd

    return nothing
end

function test_LoadScaler_sanity_checks()
    data = OPFGenerator.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))

    # Invalid noise type
    options = Dict()
    @test_throws ErrorException LoadScaler(data, options)  # missing key
    options["noise_type"] = "InvalidNoiseType"
    @test_throws ErrorException LoadScaler(data, options)  # key exists but bad value

    # Missing or invalid global parameters
    options["noise_type"] = "ScaledLogNormal"
    options["l"] = 0.8
    options["u"] = 1.2
    for v in [Inf, NaN, missing, nothing, 1+im]
        options["l"] = v
        @test_throws ErrorException LoadScaler(data, options)  # bad `l`
        options["l"] = 0.8

        options["u"] = v
        @test_throws ErrorException LoadScaler(data, options)  # bad `u`
        options["u"] = 1.2
    end
    options["l"] = 1.3
    @test_throws ErrorException LoadScaler(data, options)  # l > u
    options["l"] = 0.8

    # Invalid sigma values
    for σ in [Inf, -1, im, ["1", "2"], ones(2, 2), "0.05", [0.05, Inf]]
        options["sigma"] = σ
        @test_throws ErrorException LoadScaler(data, options)
    end

    return nothing
end

function test_sampler()
    data = OPFGenerator.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))
    _data = deepcopy(data)  # keep a deepcopy nearby
    sampler_config = Dict(
        "load" => Dict(
            "noise_type" => "ScaledLogNormal",
            "l" => 0.8,
            "u" => 1.2,
            "sigma" => 0.05,        
        )
    )
    
    rng = MersenneTwister(42)
    opf_sampler  = SimpleOPFSampler(data, sampler_config)
    data1 = rand(rng, opf_sampler)

    # No side-effect checks
    @test data == _data   # initial data should not have been modified
    @test data !== data1  # new data should be a different dictionary

    _test_ieee14_LogNormal_s42(data1)

    # Same RNG and seed should give the same data
    data2 = rand(MersenneTwister(42), opf_sampler)
    @test data2 == data1

    return nothing
end

function test_nminus1_sampler()
    data = OPFGenerator.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))
    sampler_config = Dict{String,Any}(
        "load" => Dict(
            "noise_type" => "ScaledLogNormal",
            "l" => 0.8,
            "u" => 1.2,
            "sigma" => 0.05,        
        ),
        "status" => Dict(
            "type"=> "NMinus1",
        )
    )

    rng = MersenneTwister(42)
    opf_sampler  = SimpleOPFSampler(data, sampler_config)

    data1 = rand(rng, opf_sampler)

    # all generators should be enabled
    expected_gen_status = ones(Bool, data.G)
    @test data1.gen_status == expected_gen_status

    # branch 1 should be disabled
    expected_br_status = ones(Bool, data.E)
    expected_br_status[1] = 0
    @test data1.branch_status == expected_br_status

    data2 = rand(MersenneTwister(42), opf_sampler)
    @test data2 == data1

    sampler_config["status"]["type"] = "error"
    @test_throws ErrorException SimpleOPFSampler(data, sampler_config)

    rng2 = MersenneTwister(4)
    data3 = rand(rng2, opf_sampler)
    
    # generator 3 should be disabled
    expected_gen_status = ones(Bool, data.G)
    expected_gen_status[3] = 0
    @test data3.gen_status == expected_gen_status

    # all branches should be enabled
    expected_br_status = ones(Bool, data.E)
    @test data3.branch_status == expected_br_status

    return nothing
end

function _test_ieee14_LogNormal_s42(data)
    data0 = OPFGenerator.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))

    # Check that the sampled data dictionary only has different loads/reserves
    # Buses, generators, etc... should not have changed
    # for (k, v) in data0
    for k in fieldnames(OPFGenerator.OPFData)
        v0 = getfield(data0, k)
        v = getfield(data, k)
        if k ∉ [
            :pd, :qd,
            :rmin, :rmax, :reserve_requirement,
            :branch_status, :gen_status
        ]
            @test v0 == v
        end
    end

    # Check sampled active / reactive power loads
    # The numerical values below were generated as follows, on 09/30/2024 on a RHEL 9.4 Linux machine:
    # * PGLib v21.07 case `14_ieee`, in basic network format
    # * The random number generator MersenneTwister(42)
    # * ScaledLogNormal load scaler with [0.8, 1.2] range and σ=0.05
    # ⚠ this test will fail if either condition is met
    #   * the initial data dictionary is changed
    #   * the underlying RNG or seed is changed
    #   * the load sampler config is changed
    _pd = [
        0.21477996272972988,
        0.9546062298568199,
        0.47654988934858145,
        0.08406264413456561,
        0.10703861785986431,
        0.29162856400218445,
        0.09179453210074041,
        0.031037135295027923,
        0.06490828337508349,
        0.14421853133555987,
        0.1522058057935015,
    ]
    _qd = [
        0.12570071551463455,
        0.1925426578267471,
        -0.038881685532624846,
        0.01769739876517171,
        0.071677645888302,
        0.1641028529639411,
        0.059156476242699374,
        0.01596195529458579,
        0.01702512350821862,
        0.061960554203425715,
        0.05107577375620856,
    ]
    @test data.pd ≈ _pd
    @test data.qd ≈ _qd

    @test data.reserve_requirement == 0.0
    @test data.rmin == zeros(length(data.rmin))
    @test data.rmax == zeros(length(data.rmax))

    @test all(data.gen_status)
    @test all(data.branch_status)
    return nothing
end

function test_inplace_sampler()
    data = OPFGenerator.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))
    sampler_config = Dict(
        "load" => Dict(
            "noise_type" => "ScaledLogNormal",
            "l" => 0.8,
            "u" => 1.2,
            "sigma" => 0.05,        
        )
    )

    rng = MersenneTwister(42)
    opf_sampler  = SimpleOPFSampler(data, sampler_config)
    rand!(rng, opf_sampler, data)

    _test_ieee14_LogNormal_s42(data)

    return nothing
end

function test_sampler_script()
    sampler_script = joinpath(@__DIR__, "..", "exp", "sampler.jl")
    temp_dir = mktempdir()
    config = Dict(
        "ref" => "pglib_opf_case14_ieee",
        "export_dir" => temp_dir,
        "sampler" => Dict(
            "load" => Dict(
                "noise_type" => "ScaledLogNormal",
                "l" => 0.6,
                "u" => 0.8,
                "sigma" => 0.05,
            ),
            "reserve" => Dict( # tiny reserve requirement
                "type" => "E2ELR",
                "l" => 0.0,
                "u" => 0.1,
                "factor" => 5.0,
            )
        ),
        "OPF" => Dict(
            "DCOPF" => Dict(
                "type" => "DCOPF",
                "solver" => Dict(
                    "name" => "Clarabel",
                )
            ),
            "ACOPF" => Dict(
                "type" => "ACOPF",
                "solver" => Dict(
                    "name" => "Ipopt",
                    "attributes" => Dict(
                        "tol" => 1e-6,
                    )
                )
            ),
            "ED" => Dict(
                "type" => "EconomicDispatch",
                "solver" => Dict(
                    "name" => "Clarabel",
                )
            ),
            "SOCOPF" => Dict(
                "type" => "SOCOPF",
                "solver" => Dict(
                    "name" => "Clarabel",
                )
            ),
            "SOCOPF128" => Dict(
                "type" => "SOCOPF",
                "solver" => Dict(
                    "name" => "Clarabel128",
                    "attributes" => Dict(
                        "max_iter" => 2000,
                        "max_step_fraction" => 0.995,
                        "equilibrate_enable" => true,
                        "tol_gap_abs" => 1e-12,
                        "tol_gap_rel" => 1e-12,
                        "tol_feas" => 1e-12,
                        "tol_infeas_rel" => 1e-12,
                        "tol_ktratio" => 1e-10,
                        "reduced_tol_gap_abs" => 1e-8,
                        "reduced_tol_gap_rel" => 1e-8,
                        "reduced_tol_feas" => 1e-8,
                        "reduced_tol_infeas_abs" => 1e-8,
                        "reduced_tol_infeas_rel" => 1e-8,
                        "reduced_tol_ktratio" => 1e-7,
                        "static_regularization_enable" => false,
                        "dynamic_regularization_enable" => true,
                        "dynamic_regularization_eps" => 1e-28,
                        "dynamic_regularization_delta" => 1e-14,
                        "iterative_refinement_reltol" => 1e-18,
                        "iterative_refinement_abstol" => 1e-18,
                    )
                )
            )
        )
    )

    config_file = joinpath(temp_dir, "config.toml")
    open(config_file, "w") do io
        TOML.print(io, config)
    end

    caseref = config["ref"]
    smin, smax = 1, 4
    proc = run(setenv(`$(joinpath(Sys.BINDIR, "julia")) --project=. $sampler_script $config_file $smin $smax`, dir=joinpath(@__DIR__, "..")))

    @test success(proc)

    OPFs = collect(keys(config["OPF"]))

    h5_dir = joinpath(@__DIR__, "..", config["export_dir"], "res_h5")

    @test isdir(h5_dir)

    input_file_path = joinpath(h5_dir, "$(caseref)_input_s$smin-s$smax.h5")
    @test isfile(input_file_path)
    # Check that input data file is structured as expected
    h5open(input_file_path, "r") do h5
        @test haskey(h5, "data")

        @test haskey(h5["data"], "pd")
        @test size(h5["data"]["pd"]) == (11, 4)
        @test eltype(h5["data"]["pd"]) == Float64
        @test haskey(h5["data"], "qd")
        @test size(h5["data"]["pd"]) == (11, 4)
        @test eltype(h5["data"]["qd"]) == Float64

        @test haskey(h5["data"], "branch_status")
        @test size(h5["data"]["branch_status"]) == (20, 4)
        @test eltype(h5["data"]["branch_status"]) == Bool
        @test haskey(h5["data"], "gen_status")
        @test size(h5["data"]["gen_status"]) == (5, 4)
        @test eltype(h5["data"]["gen_status"]) == Bool

        @test haskey(h5["data"], "reserve_requirement")
        @test size(h5["data"]["reserve_requirement"]) == (4,)
        @test eltype(h5["data"]["reserve_requirement"]) == Float64
        @test haskey(h5["data"], "rmin")
        @test size(h5["data"]["rmin"]) == (5,4)
        @test eltype(h5["data"]["rmin"]) == Float64
        @test haskey(h5["data"], "rmax")
        @test size(h5["data"]["rmax"]) == (5,4)
        @test eltype(h5["data"]["rmax"]) == Float64
    end

    h5_paths = [
        joinpath(h5_dir, "$(caseref)_$(opf)_s$smin-s$smax.h5")
        for opf in OPFs
    ]
    @test all(isfile.(h5_paths))
    
    for h5_path in h5_paths
        h5open(h5_path, "r") do h5
            @test haskey(h5, "meta")
            @test haskey(h5, "primal")
            @test haskey(h5, "dual")
            n_seed = length(h5["meta"]["seed"])
            @test n_seed == smax - smin + 1
            for i in 1:n_seed
                @test h5["meta"]["termination_status"][i] ∈ ["OPTIMAL", "ALMOST_OPTIMAL", "LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED"]
                @test h5["meta"]["primal_status"][i] ∈ ["FEASIBLE_POINT", "NEARLY_FEASIBLE_POINT"]
                @test h5["meta"]["dual_status"][i] ∈ ["FEASIBLE_POINT", "NEARLY_FEASIBLE_POINT"]
                @test h5["meta"]["seed"][i] == smin + i - 1
            end
        end
    end
end

@testset "Sampler" begin
    @testset test_glocal()
    @testset test_ScaledLogNormal()
    @testset test_ScaledUniform()
    @testset test_LoadScaler()
    @testset test_LoadScaler_sanity_checks()
    @testset test_sampler()
    @testset test_nminus1_sampler()
    @testset test_inplace_sampler()
    @testset test_sampler_script()
end