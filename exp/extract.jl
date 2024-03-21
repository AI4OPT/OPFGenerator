using Base.Iterators
using Base.Threads
using TOML

using HDF5
using ProgressMeter

using PowerModels
PowerModels.silence()
using OPFGenerator

if abspath(PROGRAM_FILE) == @__FILE__
    configfile = ARGS[1]
    config = TOML.parsefile(configfile)
    if get(config, "no_json", false)
        @info "Skipping JSON to HDF5 conversion since `no_json` is `true` in the config."
        exit()
    end
    batch_size = parse(Int, get(ARGS, 2, "1024"))
    OPFGenerator.parse_jsons(config;
        show_progress=false,
        batch_size,
        force_process=false,
    )
end
