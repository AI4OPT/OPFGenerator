# OPFGenerator SLURM Job Submission and Output Guide

This guide provides instructions on how to use the `submit_jobs.jl` script to submit jobs to a SLURM HPC cluster for generating OPF problem instances with OPFGenerator, as well as details on the output structure.

## Configuration

Create a TOML configuration file with the following options:

| Option | Description | Required | Default |
| ------ | ----------- | -------- | ------- |
| ref | PGLib case name | Yes | - |
| export_dir | Directory for results and intermediate files | Yes | - |
| `slurm.n_samples` | Number of samples to generate | Yes | - |
| `slurm.n_jobs` | Number of jobs to submit | Yes | - |
| `slurm.queue` | SLURM queue to submit jobs to | Yes | - |
| `slurm.charge_account` | Charge account for job submission | Yes | - |
| `slurm.minibatch_size` | Size of each minibatch | No | 16 |
| `slurm.sampler_memory` | Memory per CPU for sampler job | No | "8gb" |
| `slurm.ref_memory` | Total memory for reference job | No | `sampler_memory` |
| `slurm.sysimage_memory` | Total memory for sysimage job | No | "16gb" |
| `slurm.extract_memory` | Total memory for extract job | No | "16gb" |
| `slurm.julia_bin` | Julia command to use | No | `julia --sysimage=app/julia.so` |

## Usage

1. Create a TOML configuration file with the required options.
2. Run the `submit_jobs.jl` script with the path to the configuration file:
   ```bash
   julia --project=. slurm/submit_jobs.jl path/to/config.toml
   ```
3. Follow the printed instructions to submit the jobs to the SLURM queue.

## Output Structure

The script will organize the output into the following directories within the specified `export_dir`:

- `slurm`: Contains generated SLURM job files and submission scripts created by the `submit_jobs.jl` script.
- `res_json`: Stores JSON files with individual instance and solution data created during the sampler job.
- `res_h5`: Stores semi-aggregated HDF5 files created during the extract job.

The script will generate the following job files in the `slurm` directory:

- `sysimage.sbatch`: Job file for creating a Julia system image.
- `ref.sbatch`: Job file for solving the reference instance.
- `sampler.sbatch`: Job file for solving the instances.
- `extract.sbatch`: Job file for compiling the final HDF5 file.
- `submit.sh`: Script to submit all the above jobs to the SLURM queue.

After generating the job files, the script will output a command to run that submits the jobs. Copy+paste and run this command to start the dataset generation process. The script uses dependencies to run the jobs in the correct order -> `sysimage` -> `ref` -> `sampler` -> `extract`.

## Cleanup Script

After the job completion, you can use the `cleanup.jl` script to delete the intermediate files. This script will prompt you to confirm the deletion of the `slurm`, `res_json`, and `res_h5` directories. To run the cleanup script, execute:
```bash
julia --project=. slurm/cleanup.jl path/to/config.toml
```
You will be asked to confirm the deletion of each directory individually. This will not delete the files in the `export_dir` directory (the final results).

## Notes

- The script will check for missing instances and only generate jobs for those that are needed. Thus it is safe to re-run the script (before or after running `cleanup.jl`) if some jobs fail.
- The default memory settings won't be enough for large cases. Adjust these settings in the configuration file if necessary.