using Base.Iterators
using Base.Threads
using ProgressMeter
using HDF5

using PowerModels
PowerModels.silence()
using PGLib
using ACOPFGenerator

if abspath(PROGRAM_FILE) == @__FILE__
    result_folder = ARGS[1]
    export_folder = ARGS[2]
    batch_size    = parse(Int, ARGS[3])
    ACOPFGenerator.parse_result_folder(result_folder, export_folder;
        show_progress=false,
        batch_size,
        force_process=false,
    )
end
