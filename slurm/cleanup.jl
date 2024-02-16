using Base.Iterators
using Base.Threads
using TOML

config_file = ARGS[1]
config = TOML.parsefile(config_file)

case = config["ref"]
result_dir = config["export_dir"]

datasetname = splitdir(result_dir)[end]

opfgenerator_dir = "$(@__DIR__)/../"
name_dir = joinpath(opfgenerator_dir, "exp", datasetname)


print("Delete $name_dir? This contains job files, logs, etc. [Y/[n]] ")
if readline()[1] == 'Y'
    @info("Deleting $name_dir")
    rm(name_dir; recursive=true)
else
    @info("Not deleting $name_dir - input was not 'Y'")
end

print("Delete $(joinpath(result_dir, "res_json"))? This contains individual result files [Y/[n]] ")
if readline()[1] == 'Y'
    @info("Deleting $(joinpath(result_dir, "res_json"))")
    rm(joinpath(result_dir, "res_json"); recursive=true)
else
    @info("Not deleting $(joinpath(result_dir, "res_json")) - input was not 'Y'")
end

print("Delete $(joinpath(result_dir, "res_h5"))? This contains semi-aggregated result files [Y/[n]] ")
if readline()[1] == 'Y'
    @info("Deleting $(joinpath(result_dir, "res_h5"))")
    rm(joinpath(result_dir, "res_h5"); recursive=true)
else
    @info("Not deleting $(joinpath(result_dir, "res_h5")) - input was not 'Y'")
end
