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

function test_LoadScaler_sanity_checks()
    data = make_basic_network(pglib("pglib_opf_case14_ieee"))

    # Test potential issues
    # Not basic data
    data["basic_network"] = false
    @test_throws ErrorException LoadScaler(data, Dict())
    data["basic_network"] = true

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

function test_sampler()
    data = make_basic_network(pglib("pglib_opf_case14_ieee"))
    _data = deepcopy(data)  # keep a deepcopy nearby
    sampler_config = Dict(
        "load" => Dict(
            "noise_type" => "ScaledLogNormal",
            "l" => 0.8,
            "u" => 1.2,
            "sigma" => 0.05,        
        )
    )
    
    rng = StableRNG(42)
    opf_sampler  = SimpleOPFSampler(data, sampler_config)
    data1 = rand(rng, opf_sampler)

    # No side-effect checks
    @test data == _data   # initial data should not have been modified
    @test data !== data1  # new data should be a different dictionary

    _test_ieee14_LogNormal_s42(data1)

    # Same RNG and seed should give the same data
    data2 = rand(StableRNG(42), opf_sampler)
    @test data2 == data1

    return nothing
end

function _test_ieee14_LogNormal_s42(data)
    data0 = make_basic_network(pglib("pglib_opf_case14_ieee"))

    # Check that the sampled data dictionary only has different loads/reserves
    # Buses, generators, etc... should not have changed
    for (k, v) in data0
        if k == "gen"
            @test all(data[k][i][kk] == v[i][kk] for (i, gen) in v for (kk, vv) in gen if kk ∉ ["rmin", "rmax"])
        elseif k == "load"
            @test all(data[k][i][kk] == v[i][kk] for (i, load) in v for (kk, vv) in load if kk ∉ ["pd", "qd"])
        else
            @test data[k] == v
        end
    end

    # Check sampled active / reactive power loads
    # The numerical values below were generated as follows, on 03/06/2024 on a Linux machine:
    # * PGLib v21.07 case `14_ieee`, in basic network format
    # * The random number generator StableRNG(42)
    # * ScaledLogNormal load scaler with [0.8, 1.2] range and σ=0.05
    # ⚠ this test will fail if either condition is met
    #   * the initial data dictionary is changed
    #   * the underlying RNG or seed is changed
    #   * the load sampler config is changed
    _pd = [
        0.22876645866010775,
        1.0401653977680627,
        0.5261234671438964,
        0.07884509050225688,
        0.11947946005047866,
        0.28901062150438117,
        0.08917335057505481,
        0.03943021890027896,
        0.06710035789805892,
        0.12819215399981734,
        0.15983091550108355,
    ]
    _qd = [
        0.1338863605983119,
        0.20979981483644575,
       -0.042926391670736315,
        0.016598966421527767,
        0.08000856699808839,
        0.1626297056600925,
        0.05746727037059088,
        0.020278398291572037,
        0.0176000938749007,
        0.05507514764436597,
        0.05363453540304818,
    ]
    for (i, (p, q)) in enumerate(zip(_pd, _qd))
        @test data["load"]["$i"]["pd"] ≈ p
        @test data["load"]["$i"]["qd"] ≈ q
    end

    # all reserves should be zero
    for i in 1:length(data["gen"])
        @test data["gen"]["$i"]["rmin"] == 0.0
        @test data["gen"]["$i"]["rmax"] == 0.0
    end

    return nothing
end

function test_inplace_sampler()
    data = make_basic_network(pglib("pglib_opf_case14_ieee"))
    sampler_config = Dict(
        "load" => Dict(
            "noise_type" => "ScaledLogNormal",
            "l" => 0.8,
            "u" => 1.2,
            "sigma" => 0.05,        
        )
    )

    rng = StableRNG(42)
    opf_sampler  = SimpleOPFSampler(data, sampler_config)
    rand!(rng, opf_sampler, data)

    _test_ieee14_LogNormal_s42(data)

    return nothing
end

function test_update()
    data1 = make_basic_network(pglib("pglib_opf_case14_ieee"))
    sampler_config = Dict(
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
    )

    rng = StableRNG(42)
    opf_sampler  = SimpleOPFSampler(data1, sampler_config)
    data2 = rand(rng, opf_sampler)

    for OPF in OPFGenerator.SUPPORTED_OPF_MODELS
        solver = OPT_SOLVERS[OPF]

        opf1 = OPFGenerator.build_opf(OPF, data1, solver)
        OPFGenerator.solve!(opf1)
        OPFGenerator.update!(opf1, data2)
        
        opf2 = OPFGenerator.build_opf(OPF, data2, solver)
        
        _test_update(OPF, opf1, opf2)

        OPFGenerator.solve!(opf1)
        res1 = OPFGenerator.extract_result(opf1)

        OPFGenerator.solve!(opf2)
        res2 = OPFGenerator.extract_result(opf2)

        _test_update(OPF, opf1, opf2)

        @test res1["termination_status"] ∈ [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
        @test res2["termination_status"] ∈ [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
        @test res1["objective"] ≈ res2["objective"]
    end

    return nothing
end


function _test_update(::Type{PM.DCPPowerModel}, opf1, opf2)
    @test all(normalized_rhs.(opf1.model[:kirchhoff]) .== normalized_rhs.(opf2.model[:kirchhoff]))
end

function _test_update(::Type{OPF}, opf1, opf2) where {OPF <: Union{PM.ACPPowerModel, PM.SOCWRPowerModel, PM.SOCWRConicPowerModel}}
    @test all(normalized_rhs.(opf1.model[:kirchhoff_active]) .== normalized_rhs.(opf2.model[:kirchhoff_active]))
    @test all(normalized_rhs.(opf1.model[:kirchhoff_reactive]) .== normalized_rhs.(opf2.model[:kirchhoff_reactive]))
end

function _test_update(::Type{OPFGenerator.EconomicDispatch}, opf1, opf2)
    @test normalized_rhs(opf1.model[:power_balance]) == normalized_rhs(opf2.model[:power_balance])
    @test normalized_rhs(opf1.model[:reserve_requirement]) == normalized_rhs(opf2.model[:reserve_requirement])
    @test all(upper_bound.(opf1.model[:r]) .== upper_bound.(opf2.model[:r]))
    @test all(lower_bound.(opf1.model[:r]) .== lower_bound.(opf2.model[:r]))
    @test all(opf1.model.ext[:tracked_branches] .== opf2.model.ext[:tracked_branches])
    @test all(
        [
            normalized_rhs(opf1.model[:ptdf_flow][i]) == normalized_rhs(opf2.model[:ptdf_flow][i])
            for i in findall(opf1.model.ext[:tracked_branches])
        ]
    )
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
                "type" => "DCPPowerModel",
                "solver" => Dict(
                    "name" => "Clarabel",
                )
            ),
            "ACOPF" => Dict(
                "type" => "ACPPowerModel",
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
            "ED_noniterative" => Dict(
                "type" => "EconomicDispatch",
                "kwargs" => Dict(
                    "iterative_ptdf" => false,
                ),
                "solver" => Dict(
                    "name" => "Clarabel",
                )
            ),
            "SOCWRConic" => Dict(
                "type" => "SOCWRConicPowerModel",
                "solver" => Dict(
                    "name" => "Clarabel",
                )
            ),
            "SOCWRConic128" => Dict(
                "type" => "SOCWRConicPowerModel",
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

    @test isfile(joinpath(h5_dir, "$(caseref)_input_s$smin-s$smax.h5"))

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
                @test h5["meta"]["termination_status"][i] ∈ ["OPTIMAL", "LOCALLY_SOLVED"]
                @test h5["meta"]["primal_status"][i] == "FEASIBLE_POINT"
                @test h5["meta"]["dual_status"][i] == "FEASIBLE_POINT"
                @test h5["meta"]["seed"][i] == smin + i - 1
            end
        end
    end
end

@testset "Sampler" begin
    @testset test_glocal()
    @testset test_ScaledLogNormal()
    @testset test_ScaledUniform()
    @testset test_LoadScaler_sanity_checks()
end