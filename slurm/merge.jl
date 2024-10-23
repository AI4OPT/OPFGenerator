using OPFGenerator
using HDF5
using JSON
using TOML
using Base.Threads
using Random

main(fconfig::AbstractString) = main(TOML.parsefile(fconfig))

function main(config::Dict)
    export_dir = pop!(config, "export_dir")
    casename = config["ref"]
    all_h5_files = filter(endswith(".h5"), readdir(joinpath(export_dir, "res_h5"), join=true))

    pop!(config, "slurm")
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

        if dataset_name == "input"
            OPFGenerator.save_h5(joinpath(export_dir, "input.h5"), D)
        else
            mkpath(joinpath(export_dir, dataset_name))
            OPFGenerator.save_h5(joinpath(export_dir, dataset_name, "meta.h5"), D["meta"])
            OPFGenerator.save_h5(joinpath(export_dir, dataset_name, "primal.h5"), D["primal"])
            OPFGenerator.save_h5(joinpath(export_dir, dataset_name, "dual.h5"), D["dual"])
        end

        GC.gc()
    end

    rm(joinpath(export_dir, "res_h5"), recursive=true)

    rng = MersenneTwister(42)

    feasible = nothing
    n_samples = nothing

    # find inputs that are feasible for all OPFs
    for opf in OPFs
        h5open(joinpath(export_dir, opf, "meta.h5"), "r") do f
            if isnothing(feasible)
                n_samples = length(read(f["termination_status"]))
                feasible = trues(n_samples)
            end
            feasible .&= (String.(read(f["termination_status"])) .== "OPTIMAL") .| (String.(read(f["termination_status"])) .== "LOCALLY_SOLVED")
            feasible .&= (String.(read(f["primal_status"])) .== "FEASIBLE_POINT")
            feasible .&= (String.(read(f["dual_status"])) .== "FEASIBLE_POINT")
        end
    end

    n_feasible = sum(feasible)
    @info "Number of feasible samples: $n_feasible / $n_samples"

    feasible_idx = findall(feasible)
    infeasible_idx = sort(findall(.!feasible))

    shuffle!(rng, feasible_idx)

    n_train = Int(floor(0.8 * n_feasible))
    train_idx = sort(feasible_idx[1:n_train])
    test_idx = sort(feasible_idx[n_train+1:end])

    @info "Number of train samples: $(length(train_idx))"
    @info "Number of test samples: $(length(test_idx))"

    mkpath(joinpath(export_dir, "train"))
    mkpath(joinpath(export_dir, "test"))
    mkpath(joinpath(export_dir, "infeasible"))

    for opf in OPFs
        for file in ["primal", "dual", "meta"]
            h5open(joinpath(export_dir, opf, "$file.h5"), "r") do f
                for (idx, split) in zip([train_idx, test_idx, infeasible_idx], ["train", "test", "infeasible"])
                    mkpath(joinpath(export_dir, split, opf))
                    h5open(joinpath(export_dir, split, opf, "$file.h5"), "w") do g
                        for k in keys(f)
                            val = read(f[k])
                            g[k] = if k != "config" collect(selectdim(val, ndims(val), idx)) else val end
                        end
                    end
                end
            end
            rm(joinpath(export_dir, opf, "$file.h5"))
        end
        rm(joinpath(export_dir, opf))
    end
    # same for input
    h5open(joinpath(export_dir, "input.h5"), "r") do f
        for (idx, split) in zip([train_idx, test_idx, infeasible_idx], ["train", "test", "infeasible"])
            h5open(joinpath(export_dir, split, "input.h5"), "w") do g
                create_group(g, "data")
                for k in keys(f["data"])
                    val = read(f["data"][k])
                    g["data"][k] = collect(selectdim(val, ndims(val), idx))
                end

                create_group(g, "meta")
                for k in keys(f["meta"])
                    val = read(f["meta"][k])
                    g["meta"][k] = if k != "config" collect(selectdim(val, ndims(val), idx)) else val end
                end
            end
        end
    end
    rm(joinpath(export_dir, "input.h5"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    fconfig = ARGS[1]
    # TODO: check ARGS and that fconfig is valid path
    main(fconfig)
    exit(0)
end