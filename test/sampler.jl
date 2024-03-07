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

# function check_results_file(seed, filepath)
#     seed_status_dict = Dict(
#         1 => Dict(
#             "DCOPF" => "OPTIMAL",
#             "SOCWRConic" => "OPTIMAL",
#             "ACOPF" => "LOCALLY_SOLVED",
#         ),
#         2 => Dict(
#             "DCOPF" => "INFEASIBLE",
#             "SOCWRConic" => "INFEASIBLE",
#             "ACOPF" => "LOCALLY_INFEASIBLE",
#         ),
#         3 => Dict(
#             "DCOPF" => "OPTIMAL",
#             "SOCWRConic" => "OPTIMAL",
#             "ACOPF" => "LOCALLY_SOLVED",
#         ),
#         4 => Dict(
#             "DCOPF" => "OPTIMAL",
#             "SOCWRConic" => "INFEASIBLE",
#             "ACOPF" => "LOCALLY_INFEASIBLE",
#         ),
#         5 => Dict(
#             "DCOPF" => "OPTIMAL",
#             "SOCWRConic" => "OPTIMAL",
#             "ACOPF" => "LOCALLY_SOLVED",
#         ),
#         6 => Dict(
#             "DCOPF" => "OPTIMAL",
#             "SOCWRConic" => "SLOW_PROGRESS",
#             "ACOPF" => "LOCALLY_INFEASIBLE",
#         ),
#         7 => Dict(
#             "DCOPF" => "INFEASIBLE",
#             "SOCWRConic" => "SLOW_PROGRESS",
#             "ACOPF" => "LOCALLY_INFEASIBLE",
#         ),
#         8 => Dict(
#             "DCOPF" => "OPTIMAL",
#             "SOCWRConic" => "SLOW_PROGRESS",
#             "ACOPF" => "LOCALLY_INFEASIBLE",
#         ),
#         9 => Dict(
#             "DCOPF" => "INFEASIBLE",
#             "SOCWRConic" => "SLOW_PROGRESS",
#             "ACOPF" => "LOCALLY_INFEASIBLE",
#         ),
#         10 => Dict(
#             "DCOPF" => "OPTIMAL",
#             "SOCWRConic" => "OPTIMAL",
#             "ACOPF" => "LOCALLY_SOLVED",
#         ),
#     )

#     res = load_json(filepath)
#     for (opf, status) in seed_status_dict[seed]
#         @test res[opf]["termination_status"] == status
#         # TODO: add some solve time checks
#     end
# end

function test_sampler_script()
    sampler_script = joinpath(@__DIR__, "..", "exp", "sampler.jl")
    config_file = joinpath(@__DIR__, "config.toml")
    config = TOML.parsefile(config_file)

    caseref = config["ref"]

    proc = run(`julia --project=$(joinpath(@__DIR__, "..")) $sampler_script $config_file 1 10`)

    @test success(proc)

    # test that the output files are created and termination status is as expected
    for s in 1:10
        resfile = joinpath(config["export_dir"], "res_json", "$(caseref)_s$s.json.gz")
        @test isfile(resfile)

        # check_results_file(s, resfile)
    end
end

@testset "Sampler" begin
    test_sampler()
    test_inplace_sampler()
    test_sampler_script()
end