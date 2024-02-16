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
slurm_dir = joinpath(result_dir, "slurm")

mkpath(result_dir)
mkpath(joinpath(result_dir, "res_json"))
mkpath(joinpath(result_dir, "res_h5"))
mkpath(joinpath(slurm_dir))
mkpath(joinpath(slurm_dir, "jobs"))
mkpath(joinpath(slurm_dir, "logs"))
mkpath(joinpath(slurm_dir, "slurm"))

julia_bin = "julia --sysimage=app/julia.so"

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
    open(joinpath(slurm_dir, "jobs", "jobs_$j.txt"), "w") do io
        for minibatch in partition(seed_range, b)
            smin, smax = extrema(minibatch)
            println(io, "$(julia_bin) --project=. -t1 $(sampler_script) $(config_file) $(smin) $(smax) > $(slurm_dir)/logs/$(case)_$(smin)-$(smax).log 2>&1")
        end
    end
end

sysimage_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/sysimage.sbatch"),
    (
        charge_account=charge_account,
        slurm_dir=slurm_dir,
        queue=queue,
        opfgenerator_dir=opfgenerator_dir,
        sampler_script=sampler_script,
        config_file=config_file
    )
)
open("$(slurm_dir)/sysimage.sbatch", "w") do io
    println(io, sysimage_sbatch)
end

ref_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/ref.sbatch"),
    (
        charge_account=charge_account,
        slurm_dir=slurm_dir,
        queue=queue,
        opfgenerator_dir=opfgenerator_dir,
        config_file=config_file
    )
)
open("$(slurm_dir)/ref.sbatch", "w") do io
    println(io, ref_sbatch)
end

extract_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/extract.sbatch"),
    (
        charge_account=charge_account, 
        slurm_dir=slurm_dir, 
        queue=queue, 
        opfgenerator_dir=opfgenerator_dir, 
        extract_script=extract_script, 
        config_file=config_file
    )
)
open("$(slurm_dir)/extract.sbatch", "w") do io
    println(io, extract_sbatch)
end

sampler_sbatch = Mustache.render(
    Mustache.load("$(@__DIR__)/template/sampler.sbatch"),
    (
        charge_account=charge_account,
        slurm_dir=slurm_dir,
        queue=queue,
        J=J,
        opfgenerator_dir=opfgenerator_dir
    )
)
open("$(slurm_dir)/sampler.sbatch", "w") do io
    println(io, sampler_sbatch)
end

submit_sh = Mustache.render(
    Mustache.load("$(@__DIR__)/template/submit.sh"),
    (
        exp_dir=slurm_dir
    )
)
open("$(exp_dir)/submit.sh", "w") do io
    println(io, submit_sh)
end

@info(
"""The job files have been written to $(exp_dir).
You can run the following script to submit the jobs:"""
)

println("bash $(exp_dir)/submit.sh")

@info("Check the queue from time to time (using squeue --me)
to see if any steps failed (typically due to a lack of resources).
If so, edit the sbatch files in $(exp_dir) and re-run submit.sh.")