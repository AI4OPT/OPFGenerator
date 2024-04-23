using TOML

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

    # Check that the sampled data dictionary only has different loads
    # Buses, generators, etc... should not have changed
    for (k, v) in data0
        k == "load" && continue
        @test data[k] == v
    end
    # Only active/reactive power loads should have been modified
    for (i, load) in data0["load"]
        for (k, v) in load
            (k == "pd" || k == "qd") && continue
            @test load[k] == data["load"][i][k]
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
            "l" => 0.8,
            "u" => 1.2,
            "sigma" => 0.05,        
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
        
        @test _test_update(OPF, opf1, opf2)

        OPFGenerator.solve!(opf1)
        res1 = OPFGenerator.extract_result(opf1)

        OPFGenerator.solve!(opf2)
        res2 = OPFGenerator.extract_result(opf2)

        @test res1["termination_status"] ∈ [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
        @test res2["termination_status"] ∈ [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
        @test res1["objective"] ≈ res2["objective"]
    end

    return nothing
end


function _test_update(::Type{PM.DCPPowerModel}, opf1, opf2)
    return all(normalized_rhs.(opf1.model[:kirchhoff]) .== normalized_rhs.(opf2.model[:kirchhoff]))
end

function _test_update(::Type{OPF}, opf1, opf2) where {OPF <: Union{PM.ACPPowerModel, PM.SOCWRPowerModel, PM.SOCWRConicPowerModel}}
    return all(normalized_rhs.(opf1.model[:kirchhoff_active]) .== normalized_rhs.(opf2.model[:kirchhoff_active])) &&
            all(normalized_rhs.(opf1.model[:kirchhoff_reactive]) .== normalized_rhs.(opf2.model[:kirchhoff_reactive]))
end


function test_sampler_script()
    sampler_script = joinpath(@__DIR__, "..", "exp", "sampler.jl")
    config_file = joinpath(@__DIR__, "config.toml")
    config = TOML.parsefile(config_file)

    caseref = config["ref"]
    smin, smax = 1, 4
    proc = run(setenv(`$(joinpath(Sys.BINDIR, "julia")) --project=. $sampler_script $config_file $smin $smax`, dir=joinpath(@__DIR__, "..")))

    @test success(proc)

    OPFs = collect(keys(config["OPF"]))

    h5_dir = joinpath(@__DIR__, "..", config["export_dir"], "res_h5")

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
    test_sampler()
    test_inplace_sampler()
    test_sampler_script()
    test_update()
end