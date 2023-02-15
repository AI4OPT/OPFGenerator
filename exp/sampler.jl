using Random
using StableRNGs
using TOML

using PowerModels
PowerModels.silence()
using PGLib

using ACOPFGenerator
using JuMP
using Ipopt
using MathOptSymbolicAD

function main(rng, opf_sampler, config)

    d = Dict{String,Any}()
    d["meta"] = deepcopy(config)

    # Sample one ACOPF instance...
    d["data"] = data_ = rand(rng, opf_sampler)
    # ... solve it...
    get(config["solver"], "name", "Ipopt") == "Ipopt" || error("Only Ipopt is supported as ACOPF solver.")
    solver = optimizer_with_attributes(Ipopt.Optimizer,
        "tol" => get(config["solver"], "tol", 1e-6),
        "max_wall_time" => get(config["solver"], "max_wall_time", 3600),
    )
    acopf = ACOPFGenerator.build_acopf(data_, solver)
    # Symbolic AD is most useful for large systems
    optimize!(acopf; _differentiation_backend = MathOptSymbolicAD.DefaultBackend())
    # ... and export solution
    res = ACOPFGenerator._extract_solution(acopf, data_)
    d["res"] = res

    # Done
    return d
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Load config file
    config = TOML.parsefile(ARGS[1])
    # Parse seed range from CL arguments
    smin = parse(Int, ARGS[2])
    smax = parse(Int, ARGS[3])

    export_dir = config["export_dir"]

    # Dummy run (for pre-compilation)
    data0 = make_basic_network(pglib("14_ieee"))
    opf_sampler0 = ACOPFGenerator.SimpleOPFSampler(data0, config["sampler"])
    main(StableRNG(1), opf_sampler0, config)

    # Load reference data and setup OPF sampler
    data = make_basic_network(pglib(config["ref"]))
    opf_sampler = ACOPFGenerator.SimpleOPFSampler(data, config["sampler"])
    
    # Data generation
    @info "Generating ACOPF instances for case $ref\nSeed range: [$smin, $smax]"
    for s in smin:smax
        rng = StableRNG(s)
        tgen = @elapsed d = main(rng, opf_sampler, config)
        d["meta"]["seed"] = s
        twrite = @elapsed save_json(joinpath(export_dir, config["ref"] * "_s$s.json.gz"), d)
        @info "Seed $s" tgen twrite
    end

    return nothing
end
