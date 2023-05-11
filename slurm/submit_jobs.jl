using Base.Iterators
using Base.Threads

S = 65536  # total number of instances
J = 48     # number of jobs
b = 16     # minibatch size ()
B, r = divrem(S, J)
B = B + (r > 0)

using TOML
config_file = ARGS[1]
config = TOML.parsefile(config_file)

# case = "$(case)"
# result_dir = "data/scratch/$(case)/res_json"
case = config["ref"]
result_dir = config["export_dir"]

mkpath(result_dir)
mkpath("$(@__DIR__)/../exp/$case")
mkpath("$(@__DIR__)/../exp/$case/jobs")
mkpath("$(@__DIR__)/../exp/$case/logs")
mkpath("$(@__DIR__)/../exp/$case/slurm")

julia_bin = "julia --sysimage=app/julia.so"

for (j, seed_range) in enumerate(partition(1:S, B))
    # we run from 1 + B*(j-1) to B*j
    open("$(@__DIR__)/../exp/$(case)/jobs/jobs_$j.txt", "w") do io
        for minibatch in partition(seed_range, b)
            smin, smax = extrema(minibatch)
            println(io, "$(julia_bin) --project=. -t1 $(@__DIR__)/../exp/sampler.jl $(config_file) $(smin) $(smax) > $(@__DIR__)/../exp/$(case)/logs/$(case)_$(smin)-$(smax).log 2>&1")
        end
    end
end

# identify all missing jobs
h = falses(S)
@threads for s in 1:S
    h[s] = isfile(joinpath(result_dir, "$(case)_s$(s).json.gz"))
end
@info "Missing $(S - sum(h)) instances"
missing_seeds = (1:S)[.!h]
S_ = length(missing_seeds)
B, r = divrem(S_, J)
B = B + (r > 0)
for (j, seed_range) in enumerate(partition(missing_seeds, B))
    # Update job files
    open("$(@__DIR__)/../exp/$(case)/jobs/jobs_$j.txt", "w") do io
        for minibatch in partition(seed_range, b)
            smin, smax = extrema(minibatch)
            println(io, "$(julia_bin) --project=. -t1 $(@__DIR__)/../exp/sampler.jl $(config_file) $(smin) $(smax) > $(@__DIR__)/../exp/$(case)/logs/$(case)_$(smin)-$(smax).log 2>&1")
        end
    end
end



sampler_sbatch = "#!/bin/bash
#SBATCH --job-name=OPF            # Job name
#SBATCH --account=gts-phentenryck3-ai4opt   # charge account
#SBATCH --nodes=1                           # Use one node
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=24                  # Give 24 CPUs to each task
#SBATCH --mem-per-cpu=7gb                   # Memory per processor
#SBATCH --time=08:00:00                     # Time limit hrs:min:sec
#SBATCH -o $(@__DIR__)/../exp/$(case)/slurm/OPF.%A-%a.out         # Combined output and error messages file
#SBATCH -qembers
#SBATCH --array=1-48

module load parallel
cd $(@__DIR__)/../

srun parallel -j24 --resume-failed --joblog $(@__DIR__)/../exp/$(case)/jobs/jobs_\${SLURM_ARRAY_TASK_ID}.log {} < $(@__DIR__)/../exp/$(case)/jobs/jobs_\${SLURM_ARRAY_TASK_ID}.txt"

open("$(@__DIR__)/../exp/$(case)/sampler.sbatch", "w") do io
    println(io, sampler_sbatch)
end

extract_sbatch = "#!/bin/bash
#SBATCH --job-name=extract_OPF            # Job name
#SBATCH --account=gts-phentenryck3-ai4opt   # charge account
#SBATCH --nodes=1                           # Use one node
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=24                  # Give 128 CPUs to each task
#SBATCH --mem=0                             # Memory per processor
#SBATCH --time=08:00:00                     # Time limit hrs:min:sec
#SBATCH -o $(@__DIR__)/../exp/$(case)/slurm/extract.out         # Combined output and error messages file
#SBATCH -qembers

cd $(@__DIR__)/../

julia --project=. -t24 $(@__DIR__)/../exp/extract.jl $(result_dir) $(result_dir)/../res_h5 256"

open("$(@__DIR__)/../exp/$(case)/extract.sbatch", "w") do io
    println(io, extract_sbatch)
end

ref_sbatch = "#!/bin/bash
#SBATCH --job-name=ref_OPF            # Job name
#SBATCH --account=gts-phentenryck3-ai4opt   # charge account
#SBATCH --nodes=1                           # Use one node
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=1                   # Give 24 CPUs to each task
#SBATCH --mem-per-cpu=7gb                   # Memory per processor
#SBATCH --time=01:00:00                     # Time limit hrs:min:sec
#SBATCH -o $(@__DIR__)/../exp/$(case)/slurm/ref.out         # Combined output and error messages file
#SBATCH -qembers

module load parallel
cd $(@__DIR__)/../

julia --project=. -t1 $(@__DIR__)/../slurm/make_ref.jl $(config_file)"

open("$(@__DIR__)/../exp/$(case)/ref.sbatch", "w") do io
    println(io, ref_sbatch)
end

sysimage_sbatch = "#!/bin/bash
#SBATCH --job-name=sysimage_OPF             # Job name
#SBATCH --account=gts-phentenryck3-ai4opt   # charge account
#SBATCH --nodes=1                           # Use one node
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=1                   # Give 1 CPU
#SBATCH --mem-per-cpu=7gb                   # Memory per processor
#SBATCH --time=01:00:00                     # Time limit hrs:min:sec
#SBATCH -o $(@__DIR__)/../exp/$(case)/slurm/sysimage.out         # Combined output and error messages file
#SBATCH -qembers

cd $(@__DIR__)/../

julia --project=. -t1 --trace-compile=app/precompile.jl exp/sampler.jl exp/ieee300.toml 1 1
julia --project=. slurm/make_sysimage.jl"

open("$(@__DIR__)/../exp/$(case)/sysimage.sbatch", "w") do io
    println(io, sysimage_sbatch)
end

merge_sbatch = "#!/bin/bash
#SBATCH --job-name=merge_OPF                # Job name
#SBATCH --account=gts-phentenryck3-ai4opt   # charge account
#SBATCH --nodes=1                           # Use one node
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=1                   # Give 1 CPU
#SBATCH --mem-per-cpu=7gb                   # Memory per processor #TODO: may need to increase
#SBATCH --time=01:00:00                     # Time limit hrs:min:sec
#SBATCH -o $(@__DIR__)/../exp/$(case)/slurm/merge.out         # Combined output and error messages file
#SBATCH -qembers

cd $(@__DIR__)/../

julia --project=. slurm/merge.jl $(config_file)"

sysimage_id = read(
    `sbatch $(@__DIR__)/../exp/$(case)/sysimage.sbatch`
)
ref_id = read(
    `sbatch --dependency=afterok:$sysimage_id $(@__DIR__)/../exp/$(case)/ref.sbatch`
)
solve_id = read(
    `sbatch --dependency=afterok:$ref_id $(@__DIR__)/../exp/$(case)/sampler.sbatch`
)
extract_id = read(
    `sbatch --dependency=afterok:$solve_id $(@__DIR__)/../exp/$(case)/extract.sbatch`
)