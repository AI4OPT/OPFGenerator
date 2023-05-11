using OPFGenerator
using HDF5

using TOML
config_file = ARGS[1]
config = TOML.parsefile(config_file)
casename = config["ref"]
export_dir = config["export_dir"]

@time Ds = [
    h5read(f, "/") for f in readdir("$(export_dir)/../res_h5/", join=true) if endswith(f, ".h5")
]
@time D = OPFGenerator._merge_h5(Ds...);
@time OPFGenerator._sort_h5!(D);

# export to disk
@time OPFGenerator.save_h5("$(export_dir)/../$(casename).h5", D);