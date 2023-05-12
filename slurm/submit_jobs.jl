using Base.Iterators
using Base.Threads

S = 65536  # total number of instances
J = 45     # number of jobs
b = 16     # minibatch size ()
B, r = divrem(S, J)
B = B + (r > 0)

queue = "embers"

using TOML
config_file = ARGS[1]
config = TOML.parsefile(config_file)

case = config["ref"]
name = config["name"]

result_dir = config["export_dir"]
res_json_dir = joinpath(result_dir, "res_json")
res_h5_dir = joinpath(result_dir, "res_h5")

opfgenerator_dir = "$(@__DIR__)/../"
sampler_script = joinpath(opfgenerator_dir, "exp", "sampler.jl")
extract_script = joinpath(opfgenerator_dir, "exp", "extract.jl")
name_dir = joinpath(opfgenerator_dir, "exp", name)

mkpath(result_dir)
mkpath(joinpath(result_dir, "res_json"))
mkpath(joinpath(result_dir, "res_h5"))
mkpath(joinpath(opfgenerator_dir, "exp", name))
mkpath(joinpath(opfgenerator_dir, "exp", name, "jobs"))
mkpath(joinpath(opfgenerator_dir, "exp", name, "logs"))
mkpath(joinpath(opfgenerator_dir, "exp", name, "slurm"))

julia_bin = "julia --sysimage=app/julia.so"

for (j, seed_range) in enumerate(partition(1:S, B))
    # we run from 1 + B*(j-1) to B*j
    open(joinpath(opfgenerator_dir, "exp", name, "jobs", "jobs_$j.txt"), "w") do io
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
    open(joinpath(opfgenerator_dir, "exp", name, "jobs", "jobs_$j.txt"), "w") do io
        for minibatch in partition(seed_range, b)
            smin, smax = extrema(minibatch)
            println(io, "$(julia_bin) --project=. -t1 $(sampler_script) $(config_file) $(smin) $(smax) > $(name_dir)/logs/$(case)_$(smin)-$(smax).log 2>&1")
        end
    end
end

sampler_sbatch = "#!/bin/bash
#SBATCH --job-name=OPF                      # Job name
#SBATCH --account=gts-phentenryck3-ai4opt   # charge account
#SBATCH --nodes=1                           # Use one node
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=24                  # Give 24 CPUs to each task
#SBATCH --mem-per-cpu=7gb                   # Memory per processor
#SBATCH --time=08:00:00                     # Time limit hrs:min:sec
#SBATCH -o $(name_dir)/slurm/OPF.%A-%a.out         # Combined output and error messages file
#SBATCH -q$(queue)
#SBATCH --array=1-$(J)

module load parallel

cd $(opfgenerator_dir)

srun parallel -j24 --resume-failed --joblog $(name_dir)/jobs/jobs_\${SLURM_ARRAY_TASK_ID}.log {} < $(name_dir)/jobs/jobs_\${SLURM_ARRAY_TASK_ID}.txt"

open("$(name_dir)/sampler.sbatch", "w") do io
    println(io, sampler_sbatch)
end

extract_sbatch = "#!/bin/bash
#SBATCH --job-name=extract_OPF              # Job name
#SBATCH --account=gts-phentenryck3-ai4opt   # charge account
#SBATCH --nodes=1                           # Use one node
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=24                  # Give 128 CPUs to each task
#SBATCH --mem=7gb                           # Memory per processor
#SBATCH --time=08:00:00                     # Time limit hrs:min:sec
#SBATCH -o $(name_dir)/slurm/extract.out         # Combined output and error messages file
#SBATCH -q$(queue)

cd $(opfgenerator_dir)

julia --project=. -t24 $(extract_script) $(res_json_dir) $(res_h5_dir) 256
julia --project=. slurm/merge.jl $(config_file)"

open("$(name_dir)/extract.sbatch", "w") do io
    println(io, extract_sbatch)
end

ref_sbatch = "#!/bin/bash
#SBATCH --job-name=ref_OPF                  # Job name
#SBATCH --account=gts-phentenryck3-ai4opt   # charge account
#SBATCH --nodes=1                           # Use one node
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=1                   # Give 24 CPUs to each task
#SBATCH --mem-per-cpu=7gb                   # Memory per processor
#SBATCH --time=01:00:00                     # Time limit hrs:min:sec
#SBATCH -o $(name_dir)/slurm/ref.out         # Combined output and error messages file
#SBATCH -q$(queue)

cd $(opfgenerator_dir)

julia --project=. -t1 $(opfgenerator_dir)/slurm/make_ref.jl $(config_file)"

open("$(name_dir)/ref.sbatch", "w") do io
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
#SBATCH -o $(name_dir)/slurm/sysimage.out         # Combined output and error messages file
#SBATCH -q$(queue)

cd $(opfgenerator_dir)

mkdir app

julia --project=. -t1 --trace-compile=app/precompile.jl $(sampler_script) $(config_file) 1 1

julia --project=. slurm/make_sysimage.jl"

open("$(name_dir)/sysimage.sbatch", "w") do io
    println(io, sysimage_sbatch)
end


sysimage_id = split(readchomp(
    `sbatch $(name_dir)/sysimage.sbatch` # took 12 min
))[end]
println("Submitted sysimage job with id $sysimage_id")

ref_id = split(readchomp(
    `sbatch --dependency=afterok:$sysimage_id $(name_dir)/ref.sbatch` # took 2 min with 1354pegase
))[end]
println("Submitted ref job with id $ref_id")

solve_id = split(readchomp(
    `sbatch --dependency=afterok:$ref_id $(name_dir)/sampler.sbatch`
))[end]
println("Submitted solve job with id $solve_id")

# extract_id = split(readchomp( # TODO: for some reason this errors out when submitting automatically, but works when submitting manually...
#     `sbatch --dependency=afterok:$solve_id $(name_dir)/extract.sbatch`
# ))[end]
# println("Submitted extract+merge job with id $extract_id")
@info("Run the following command to submit the extract+merge job:
sbatch --dependency=afterok:$solve_id $(name_dir)/extract.sbatch")