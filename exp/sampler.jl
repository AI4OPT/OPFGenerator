using Random
using LinearAlgebra
using StableRNGs
using TOML

using PowerModels
PowerModels.silence()
using PGLib

using OPFGenerator
using JuMP
using Ipopt
using HSL

using MathOptSymbolicAD

const LIB_COINHSL = HSL.libcoinhsl

function main(rng, opf_sampler, config)

    d = Dict{String,Any}()
    d["meta"] = deepcopy(config)

    # Sample one ACOPF instance...
    d["data"] = data_ = rand(rng, opf_sampler)
    # ... solve it...
    get(config["solver"], "name", "Ipopt") == "Ipopt" || error("Only Ipopt is supported as ACOPF solver.")
    solver = optimizer_with_attributes(Ipopt.Optimizer,
        "hsllib" => LIB_COINHSL,
        "tol" => get(config["solver"], "tol", 1e-6),
        "max_wall_time" => get(config["solver"], "max_wall_time", 3600.0),
        "linear_solver" => get(config["solver"], "linear_solver", "mumps")
    )
    acopf = OPFGenerator.build_acopf(data_, solver)
    set_silent(acopf)
    # Symbolic AD is most useful for large systems
    optimize!(acopf; _differentiation_backend = MathOptSymbolicAD.DefaultBackend())
    # ... and export solution
    acopf_res = OPFGenerator._extract_acopf_solution(acopf, data_)
    d["acopf_res"] = acopf_res

    dcopf = OPFGenerator.build_dcopf(data_, solver)
    set_silent(dcopf)
    optimize!(dcopf; _differentiation_backend = MathOptSymbolicAD.DefaultBackend())
    dcopf_res = OPFGenerator._extract_dcopf_solution(dcopf, data_)
    d["dcopf_res"] = dcopf_res

    # Done
    return d
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Load config file
    config = TOML.parsefile(ARGS[1])
    # Parse seed range from CL arguments
    smin = parse(Int, ARGS[2])
    smax = parse(Int, ARGS[3])

    export_dir = joinpath(config["export_dir"], "res_json")

    # Dummy run (for pre-compilation)
    data0 = make_basic_network(pglib("14_ieee"))
    opf_sampler0 = OPFGenerator.SimpleOPFSampler(data0, config["sampler"])
    main(StableRNG(1), opf_sampler0, config)

    # Load reference data and setup OPF sampler
    data = make_basic_network(pglib(config["ref"]))
    opf_sampler = OPFGenerator.SimpleOPFSampler(data, config["sampler"])
    
    # Data generation
    @info "Generating ACOPF & DCOPF instances for case $ref\nSeed range: [$smin, $smax]"
    for s in smin:smax
        rng = StableRNG(s)
        tgen = @elapsed d = main(rng, opf_sampler, config)
        d["meta"]["seed"] = s
        twrite = @elapsed save_json(joinpath(export_dir, config["ref"] * "_s$s.json.gz"), d)
        @info "Seed $s" tgen twrite
    end

    return nothing
end
