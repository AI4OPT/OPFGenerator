#!/usr/bin/env bash

sysimage_exists="no"
prompt_message="Do you want to create a sysimage at app/julia.so?"
default_sysimage="yes"
if [ -f "app/julia.so" ]; then
    sysimage_exists="yes"
    prompt_message="Do you want to re-create the sysimage at app/julia.so?"
    default_sysimage="no"
fi

ask_sysimage() {
    while true; do
        read -p "${prompt_message} (yes/no) [${default_sysimage}]: " yn
        yn=${yn:-$default_sysimage} # Default value if input is empty
        case $(echo "$yn" | tr '[:upper:]' '[:lower:]') in
            'y' | 'yes' ) run_sysimage="yes"; break;;
            'n' | 'no' ) run_sysimage="no"; break;;
            * ) echo "Error: Unrecognized input. Please answer 'yes' or 'no'.";;
        esac
    done
}

ask_sysimage

if [ "$run_sysimage" == "yes" ]; then
    if [ "$sysimage_exists" == "yes" ]; then
        echo "Re-creating sysimage..."
    else
        echo "Creating a sysimage..."
    fi
    sysimage_id=$(sbatch {{{:slurm_dir}}}/sysimage.sbatch | awk '{print $NF}')
    echo "Submitted sysimage job with id $sysimage_id"
    dependency_sysimage="--dependency=afterok:$sysimage_id"
else
    echo "Skipping sysimage job"
    dependency_sysimage=""
fi

ref_id=$(sbatch $dependency_sysimage {{{:slurm_dir}}}/ref.sbatch | awk '{print $NF}')
echo "Submitted ref job with id $ref_id"

solve_id=$(sbatch --dependency=afterok:$ref_id {{{:slurm_dir}}}/sampler.sbatch | awk '{print $NF}')
echo "Submitted solve job with id $solve_id"

merge_id=$(sbatch --dependency=afterok:$solve_id {{{:slurm_dir}}}/extract.sbatch | awk '{print $NF}')
echo "Submitted extract+merge job with id $merge_id"