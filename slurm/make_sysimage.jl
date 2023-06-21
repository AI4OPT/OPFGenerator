using PackageCompiler

# first run `julia --project=. -t1 --trace-compile=app/precompile.jl exp/sampler.jl exp/ieee300.toml 1 1`
create_sysimage(
    [
        "OPFGenerator",
        "Distributions",
        "HSL",
        "Ipopt",
        "JSON",
        "JuMP",
        "LinearAlgebra",
        "MathOptSymbolicAD",
        "PGLib",
        "PowerModels",
        "ProgressMeter",
        "Random",
        "StableRNGs",
        "TOML",
    ];
    sysimage_path="app/julia.so",
    precompile_statements_file="app/precompile.jl"
);