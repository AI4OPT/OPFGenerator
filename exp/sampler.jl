using Random
using StableRNGs
using TOML

using PowerModels
PowerModels.silence()
using PGLib

using ACOPFGenerator
using JuMP
using Ipopt

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
    set_silent(acopf)
    optimize!(acopf)
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

    # Load reference data
    data = make_basic_network(pglib(config["ref"]))
    # Setup OPF sampler
    opf_sampler = ACOPFGenerator.SimpleOPFSampler(data, config["sampler"])

    # Dummy run (for pre-compilation)
    main(StableRNG(42), opf_sampler, config)

    # Data generation
    @info "Generating ACOPF instances for case $ref\nSeed range: [$smin, $smax]"
    for s in smin:smax
        rng = StableRNG(s)
        d = main(rng, opf_sampler, config)
        d["meta"]["seed"] = s
        save_json(joinpath(export_dir, config["ref"] * "_s$s.json.gz"), d)
    end

    return nothing
end
