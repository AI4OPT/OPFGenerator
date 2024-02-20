# Utilities for using OPFGenerator with Slurm HPC


This folder contains some scripts and templates to run OPFGenerator on a SLURM HPC cluster.

The file `submit_jobs.jl` is the main entry point and will automatically submit dependent jobs to the cluster. It is configured by default to use Pascal's AI4OPT charge account and the `embers` queue. You can change these settings by editing the first few lines of `submit_jobs.jl`.

To launch a dataset generation job, first create a TOML file with the desired parameters. Then, run the following command from the root directory of the `OPFGenerator` repo:

```bash
julia --project=. slurm/submit_jobs.jl <path/to/config/file.toml>
```

Follow the directions from the output of the above command to submit the jobs.

Once submitted, check the queue from time to time (using `squeue --me`) to see if any steps failed (typically due to a lack of resources). If so, edit the generated sbatch files (in <export_dir>/slurm) and re-run the generated `submit.sh` script (also in <export_dir>/slurm).

When the jobs are finished, you can run the below interactive script to remove any temporary files:

```bash
julia --project=. slurm/cleanup.jl <path/to/config/file.toml>
```