using OPFGenerator
using HDF5
using TOML
using Base.Threads

config_file = ARGS[1]
config = TOML.parsefile(config_file)
name = config["name"]
export_dir = config["export_dir"]

# @time Ds = [
#     h5read(f, "/") for f in readdir("$(export_dir)/res_h5/", join=true) if endswith(f, ".h5")
# ]
function read_files(export_dir)
    h5_files = [f for f in readdir("$(export_dir)/res_h5/", join=true) if endswith(f, ".h5")]
    n_files = length(h5_files)
    Ds = Vector{Dict{String,Any}}(undef, n_files)

    @threads for i in 1:n_files
        Ds[i] = h5read(h5_files[i], "/")
    end
    return Ds
end

@time Ds = read_files(export_dir);

@time D = OPFGenerator._merge_h5(Ds...);
@time OPFGenerator._sort_h5!(D);

# export to disk
@time OPFGenerator.save_h5("$(export_dir)/$(name).h5", D);