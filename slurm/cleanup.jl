using TOML

config_file = ARGS[1]
config = TOML.parsefile(config_file)

case = config["ref"]
result_dir = config["export_dir"]
slurm_dir = joinpath(result_dir, "slurm")
json_dir = joinpath(result_dir, "res_json")
h5_dir = joinpath(result_dir, "res_h5")


print("Delete $slurm_dir? This contains job files, logs, etc. [Y/[n]] ")
response = readline()
if response == "Y" || response == "y"
    @info("Deleting $slurm_dir")
    rm(slurm_dir; recursive=true)
else
    @info("Not deleting $slurm_dir - got input '$response'")
end

print("Delete $json_dir? This contains individual result files [Y/[n]] ")
response = readline()
if response == "Y" || response == "y"
    @info("Deleting $json_dir")
    rm(json_dir; recursive=true)
else
    @info("Not deleting $json_dir - got input '$response'")
end

print("Delete $h5_dir)? This contains semi-aggregated result files [Y/[n]] ")
response = readline()
if response == "Y" || response == "y"
    @info("Deleting $h5_dir")
    rm(h5_dir; recursive=true)
else
    @info("Not deleting $h5_dir - got input '$response'")
end
