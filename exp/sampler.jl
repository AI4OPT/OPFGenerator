using Random
using LinearAlgebra
using StableRNGs
using TOML

using PowerModels
PowerModels.silence()
using PGLib

using JuMP
using Ipopt
using Mosek, MosekTools

using HSL_jll
const LIB_COINHSL = HSL_jll.libhsl_path

using MathOptSymbolicAD

using OPFGenerator

const NAME2OPTIMIZER = Dict(
    "Ipopt" => Ipopt.Optimizer,
    "Mosek" => Mosek.Optimizer,
)

function main(data, config)
    d = Dict{String,Any}()
    d["meta"] = deepcopy(config)

    # Keep track of initial data file
    d["data"] = data

    # Solve all OPF formulations
    for (opf_str, opf_solver) in config["OPF"]
        OPF = OPFGenerator.OPF2TYPE[opf_str]
        
        solver = optimizer_with_attributes(NAME2OPTIMIZER[opf_solver],
            config["solver"][opf_solver]...
        )

        opf = OPFGenerator.build_opf(OPF, data, solver)
       
        # Solve OPF model
        set_silent(opf.model)
        optimize!(opf.model; _differentiation_backend = MathOptSymbolicAD.DefaultBackend())

        res = OPFGenerator.extract_result(opf)
        d[opf_str] = res
    end

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
    mkpath(export_dir)

    # Dummy run (for pre-compilation)
    data0 = make_basic_network(pglib("14_ieee"))
    opf_sampler0 = OPFGenerator.SimpleOPFSampler(data0, config["sampler"])
    rand(StableRNG(1), opf_sampler0)
    d0 = main(data0, config)

    # Load reference data and setup OPF sampler
    data = make_basic_network(pglib(config["ref"]))
    opf_sampler = OPFGenerator.SimpleOPFSampler(data, config["sampler"])
    
    # Data generation
    @info "Generating ACOPF & DCOPF instances for case $ref\nSeed range: [$smin, $smax]"
    for s in smin:smax
        rng = StableRNG(s)
        data_ = rand(rng, opf_sampler)
        tgen = @elapsed d = main(data_, config)
        d["meta"]["seed"] = s
        twrite = @elapsed save_json(joinpath(export_dir, config["ref"] * "_s$s.json.gz"), d)
        @info "Seed $s" tgen twrite
    end

    return nothing
end
