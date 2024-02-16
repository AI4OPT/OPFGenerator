#!/bin/bash
#SBATCH --job-name=OPF                        # Job name
#SBATCH --account={{:charge_account}}         # charge account
#SBATCH --nodes=1                             # Use one node
#SBATCH --ntasks=1                            # Run a single task
#SBATCH --cpus-per-task=24                    # Give 24 CPUs to each task
#SBATCH --mem-per-cpu=7gb                     # Memory per processor
#SBATCH --time=08:00:00                       # Time limit hrs:min:sec
#SBATCH -o {{{:slurm_dir}}}/slurm/OPF.%A-%a.out  # Combined output and error messages file
#SBATCH -q{{:queue}}
#SBATCH --array=1-{{:J}}

module load parallel

cd {{{:opfgenerator_dir}}}

srun parallel -j24 --resume-failed --joblog {{{:slurm_dir}}}/jobs/jobs_${SLURM_ARRAY_TASK_ID}.log {} < {{{:slurm_dir}}}/jobs/jobs_${SLURM_ARRAY_TASK_ID}.txt