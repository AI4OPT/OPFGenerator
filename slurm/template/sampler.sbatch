#!/bin/bash
#SBATCH --job-name=OPF                        # Job name
#SBATCH --account={{:charge_account}}         # charge account
#SBATCH --nodes=1                             # Use one node
#SBATCH --ntasks=1                            # Run a single task
#SBATCH --cpus-per-task={{:cpus_per_task}}    # Give CPUs to each task
#SBATCH --mem-per-cpu={{:sampler_memory}}     # Memory per processor
#SBATCH --time=08:00:00                       # Time limit hrs:min:sec
#SBATCH -o {{{:logs_dir}}}/OPF.%A-%a.out      # Combined output and error messages file
#SBATCH -q{{:queue}}
#SBATCH --array=1-{{:J}}

. {{{:env_path}}}

cd {{{:opfgenerator_dir}}}

srun parallel -j{{:cpus_per_task}} --resume-failed --joblog {{{:jobs_dir}}}/jobs_${SLURM_ARRAY_TASK_ID}.log {} < {{{:jobs_dir}}}/jobs_${SLURM_ARRAY_TASK_ID}.txt