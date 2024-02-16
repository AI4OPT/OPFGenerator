using Base.Iterators
using Base.Threads
using Mustache
using TOML

S = 65536  # total number of instances
J = 45     # number of jobs
b = 16     # minibatch size ()
B, r = divrem(S, J)
B = B + (r > 0)

queue = "embers"
charge_account = "gts-phentenryck3-ai4opt"

config_file = ARGS[1]
config = TOML.parsefile(config_file)

case = config["ref"]
result_dir = config["export_dir"]

datasetname = splitdir(result_dir)[end]

opfgenerator_dir = "$(@__DIR__)/../"
sampler_script = joinpath(opfgenerator_dir, "exp", "sampler.jl")
extract_script = joinpath(opfgenerator_dir, "exp", "extract.jl")
name_dir = joinpath(opfgenerator_dir, "exp", datasetname)

mkpath(result_dir)
mkpath(joinpath(result_dir, "res_json"))
mkpath(joinpath(result_dir, "res_h5"))
mkpath(joinpath(opfgenerator_dir, "exp", datasetname))
mkpath(joinpath(opfgenerator_dir, "exp", datasetname, "jobs"))
mkpath(joinpath(opfgenerator_dir, "exp", datasetname, "logs"))
mkpath(joinpath(opfgenerator_dir, "exp", datasetname, "slurm"))

julia_bin = "julia --sysimage=app/julia.so"

for (j, seed_range) in enumerate(partition(1:S, B))
    # we run from 1 + B*(j-1) to B*j
    open(joinpath(opfgenerator_dir, "exp", datasetname, "jobs", "jobs_$j.txt"), "w") do io
        for minibatch in partition(seed_range, b)
            smin, smax = extrema(minibatch)
            println(io, "$(julia_bin) --project=. -t1 $(sampler_script) $(config_file) $(smin) $(smax) > $(name_dir)/logs/$(case)_$(smin)-$(smax).log 2>&1")
        end
    end
end

# identify all missing jobs
h = falses(S)
@threads for s in 1:S
    h[s] = isfile(joinpath(result_dir, "res_json", "$(case)_s$(s).json.gz"))
end
@info "Missing $(S - sum(h)) instances"
missing_seeds = (1:S)[.!h]
S_ = length(missing_seeds)
B, r = divrem(S_, J)
B = B + (r > 0)
for (j, seed_range) in enumerate(partition(missing_seeds, B))
    # Update job files
    open(joinpath(opfgenerator_dir, "exp", datasetname, "jobs", "jobs_$j.txt"), "w") do io
        for minibatch in partition(seed_range, b)
            smin, smax = extrema(minibatch)
            println(io, "$(julia_bin) --project=. -t1 $(sampler_script) $(config_file) $(smin) $(smax) > $(name_dir)/logs/$(case)_$(smin)-$(smax).log 2>&1")
        end
    end
end

sysimage_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/sysimage.sbatch"),
    (
        charge_account=charge_account,
        name_dir=name_dir,
        queue=queue,
        opfgenerator_dir=opfgenerator_dir,
        sampler_script=sampler_script,
        config_file=config_file
    )
)
open("$(name_dir)/sysimage.sbatch", "w") do io
    println(io, sysimage_sbatch)
end

ref_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/ref.sbatch"),
    (
        charge_account=charge_account,
        name_dir=name_dir,
        queue=queue,
        opfgenerator_dir=opfgenerator_dir,
        config_file=config_file
    )
)
open("$(name_dir)/ref.sbatch", "w") do io
    println(io, ref_sbatch)
end

extract_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/extract.sbatch"),
    (
        charge_account=charge_account, 
        name_dir=name_dir, 
        queue=queue, 
        opfgenerator_dir=opfgenerator_dir, 
        extract_script=extract_script, 
        config_file=config_file
    )
)
open("$(name_dir)/extract.sbatch", "w") do io
    println(io, extract_sbatch)
end

sampler_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/sampler.sbatch"),
    (
        charge_account=charge_account,
        name_dir=name_dir,
        queue=queue,
        J=J,
        opfgenerator_dir=opfgenerator_dir
    )
)
open("$(name_dir)/sampler.sbatch", "w") do io
    println(io, sampler_sbatch)
end


sysimage_id = split(readchomp(
    `sbatch $(name_dir)/sysimage.sbatch`
))[end]
println("Submitted sysimage job with id $sysimage_id")

ref_id = split(readchomp(
    `sbatch --dependency=afterok:$sysimage_id $(name_dir)/ref.sbatch`
))[end]
println("Submitted ref job with id $ref_id")

solve_id = split(readchomp(
    `sbatch --dependency=afterok:$ref_id $(name_dir)/sampler.sbatch`
))[end]
println("Submitted solve job with id $solve_id")

# TODO: for some reason this errors out when submitting automatically, but works when submitting manually...
# extract_id = split(readchomp(
#     `sbatch --dependency=afterok:$solve_id $(name_dir)/extract.sbatch`
# ))[end]
# println("Submitted extract+merge job with id $extract_id")
@info("Run the following command to submit the extract+merge job:
sbatch --dependency=afterok:$solve_id $(name_dir)/extract.sbatch")