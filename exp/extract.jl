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
    batch_size = parse(Int, get(ARGS, 2, "1024"))
    OPFGenerator.parse_jsons(config;
        show_progress=false,
        batch_size,
        force_process=false,
    )
end
