#!/usr/bin/env bash

sysimage_id=$(sbatch {{{:slurm_dir}}}/sysimage.sbatch | awk '{print $NF}')
echo "Submitted sysimage job with id $sysimage_id"

ref_id=$(sbatch --dependency=afterok:$sysimage_id {{{:slurm_dir}}}/ref.sbatch | awk '{print $NF}')
echo "Submitted ref job with id $ref_id"

solve_id=$(sbatch --dependency=afterok:$ref_id {{{:slurm_dir}}}/sampler.sbatch | awk '{print $NF}')
echo "Submitted solve job with id $solve_id"

merge_id=$(sbatch --dependency=afterok:$solve_id {{{:slurm_dir}}}/extract.sbatch | awk '{print $NF}')
echo "Submitted extract+merge job with id $merge_id"