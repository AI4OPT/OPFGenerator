using PackageCompiler

# first run `julia --project=. -t1 --trace-compile=app/precompile.jl exp/sampler.jl exp/ieee300.toml 1 1`
create_sysimage(
    [
        "OPFGenerator",
        "Random",
        "StableRNGs",
        "Distributions",
        "PowerModels",
        "PGLib",
        "JuMP",
        "Graphs",
        "ProgressMeter",
        "HDF5",
        "Ipopt",
        "Mosek",
        "MosekTools",
        "HSL_jll",
        "LinearAlgebra",
        "MathOptSymbolicAD",
        "TOML",
        "CodecBzip2",
        "CodecZlib",
        "JSON"
    ];
    sysimage_path="app/julia.so",
    precompile_statements_file="app/precompile.jl"
);