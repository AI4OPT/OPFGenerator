#!/bin/bash
#SBATCH --job-name=extract_OPF                # Job name
#SBATCH --account={{:charge_account}}         # charge account
#SBATCH --nodes=1                             # Use one node
#SBATCH --ntasks=1                            # Run a single task
#SBATCH --cpus-per-task=24                    # Give 128 CPUs to each task
#SBATCH --mem=7gb                             # Memory per processor
#SBATCH --time=08:00:00                       # Time limit hrs:min:sec
#SBATCH -o {{{:slurm_dir}}}/slurm/extract.out    # Combined output and error messages file
#SBATCH -q{{:queue}}

cd {{{:opfgenerator_dir}}}

julia --project=. -t24 {{{:extract_script}}} {{{:config_file}}}
julia --project=. -t24 slurm/merge.jl {{{:config_file}}}