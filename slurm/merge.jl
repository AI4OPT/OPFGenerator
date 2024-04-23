using OPFGenerator
using HDF5
using JSON
using TOML
using Base.Threads

main(fconfig::AbstractString) = main(TOML.parsefile(fconfig))

function main(config::Dict)
    export_dir = config["export_dir"]
    casename = config["ref"]
    all_h5_files = filter(endswith(".h5"), readdir(joinpath(export_dir, "res_h5"), join=true))
    config_str = JSON.json(config)

    # Process each dataset
    OPFs = sort(collect(keys(config["OPF"])))
    for dataset_name in ["input"; OPFs]
        rgx = Regex("$(casename)_$(dataset_name)_s\\d+-s\\d+.h5")
        res_files = filter(s -> startswith(basename(s), rgx), all_h5_files)
        nfiles = length(res_files)
        Ds = Vector{Dict{String,Any}}(undef, nfiles)
        @threads for i in 1:nfiles
            Ds[i] = h5read(res_files[i], "/")
        end
        
        # Merge minibatches, sort, and export to disk
        D = OPFGenerator._merge_h5(Ds)
        OPFGenerator._dedupe_and_sort_h5!(D)

        # Save dataset to disk
        get!(D, "meta", Dict{String,Any}())
        D["meta"]["config"] = config_str
        OPFGenerator.save_h5(joinpath(export_dir, "$(casename)_$(dataset_name).h5"), D);

        GC.gc()
    end

    return nothing    
end

if abspath(PROGRAM_FILE) == @__FILE__
    fconfig = ARGS[1]
    # TODO: check ARGS and that fconfig is valid path
    main(fconfig)
    exit(0)
end
