using Base.Iterators
using Base.Threads
using Mustache
using TOML


config_file = ARGS[1]
config = TOML.parsefile(config_file)

case = config["ref"]
result_dir = config["export_dir"]
S = config["slurm"]["n_samples"]
J = config["slurm"]["n_jobs"]
queue = config["slurm"]["queue"]
charge_account = config["slurm"]["charge_account"]

b = get(config["slurm"], "minibatch_size", 16)
sampler_memory = get(config["slurm"], "sampler_memory", "8gb")
extract_memory = get(config["slurm"], "extract_memory", "64gb")
ref_memory = get(config["slurm"], "ref_memory", sampler_memory)
sysimage_memory = get(config["slurm"], "sysimage_memory", "16gb")
julia_bin = get(config["slurm"], "julia_bin", "julia --sysimage=app/julia.so")
cpus_per_task = get(config["slurm"], "cpus_per_task", 24)
env_path = get(config["slurm"], "env_path", joinpath(@__DIR__, "template", "env.sh"))

B, r = divrem(S, J)
B = B + (r > 0)

datasetname = splitdir(result_dir)[end]

opfgenerator_dir = "$(@__DIR__)/../"
sampler_script = joinpath(opfgenerator_dir, "exp", "sampler.jl")

slurm_dir = joinpath(result_dir, "slurm")
json_dir = joinpath(result_dir, "res_json")
h5_dir = joinpath(result_dir, "res_h5")

jobs_dir = joinpath(slurm_dir, "jobs")
logs_dir = joinpath(slurm_dir, "logs")

mkpath(result_dir)
mkpath(json_dir)
mkpath(h5_dir)
mkpath(slurm_dir)
mkpath(jobs_dir)
mkpath(logs_dir)

# identify all missing jobs
h = falses(S)
@threads for s in 1:S
    h[s] = isfile(joinpath(json_dir, "$(case)_s$(s).json.gz"))
end
@info "Missing $(S - sum(h)) instances"
missing_seeds = (1:S)[.!h]
S_ = length(missing_seeds)
B, r = divrem(S_, J)
B = B + (r > 0)
for (j, seed_range) in enumerate(partition(missing_seeds, B))
    # Update job files
    open(joinpath(jobs_dir, "jobs_$j.txt"), "w") do io
        for minibatch in partition(seed_range, b)
            smin, smax = extrema(minibatch)
            println(io, "$(julia_bin) --project=. -t1 $(sampler_script) $(config_file) $(smin) $(smax) > $(logs_dir)/$(case)_$(smin)-$(smax).log 2>&1")
        end
    end
end

sysimage_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/sysimage.sbatch"),
    (
        charge_account=charge_account,
        logs_dir=logs_dir,
        queue=queue,
        opfgenerator_dir=opfgenerator_dir,
        sampler_script=sampler_script,
        config_file=config_file,
        sysimage_memory=sysimage_memory,
        env_path=env_path,
    )
)
open("$(slurm_dir)/sysimage.sbatch", "w") do io
    println(io, sysimage_sbatch)
end

ref_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/ref.sbatch"),
    (
        charge_account=charge_account,
        logs_dir=logs_dir,
        queue=queue,
        opfgenerator_dir=opfgenerator_dir,
        config_file=config_file,
        ref_memory=ref_memory,
        env_path=env_path,
    )
)
open("$(slurm_dir)/ref.sbatch", "w") do io
    println(io, ref_sbatch)
end

extract_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/extract.sbatch"),
    (
        charge_account=charge_account, 
        logs_dir=logs_dir, 
        queue=queue, 
        opfgenerator_dir=opfgenerator_dir,
        config_file=config_file,
        extract_memory=extract_memory,
        cpus_per_task=cpus_per_task,
        env_path=env_path,
    )
)
open("$(slurm_dir)/extract.sbatch", "w") do io
    println(io, extract_sbatch)
end

sampler_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/sampler.sbatch"),
    (
        charge_account=charge_account,
        logs_dir=logs_dir,
        jobs_dir=jobs_dir,
        queue=queue,
        J=J,
        opfgenerator_dir=opfgenerator_dir,
        sampler_memory=sampler_memory,
        cpus_per_task=cpus_per_task,
        env_path=env_path,
    )
)
open("$(slurm_dir)/sampler.sbatch", "w") do io
    println(io, sampler_sbatch)
end

submit_sh = Mustache.render(
    Mustache.load("$(@__DIR__)/template/submit.sh"),
    (
        slurm_dir=slurm_dir
    )
)
open("$(slurm_dir)/submit.sh", "w") do io
    println(io, submit_sh)
end

run(`chmod +x $(slurm_dir)/submit.sh`)

@info(
"""The job files have been written to $(slurm_dir).
Run the below command to submit them to the queue."""
)

command = "bash $(slurm_dir)/submit.sh"
command_len = length(command)
println("↓↓↓" * repeat("↓", command_len) * "↓↓↓")
println("\n   " * command)
println("\n↑↑↑" * repeat("↑", command_len) * "↑↑↑")