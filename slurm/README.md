# Utilities for using OPFGenerator with Slurm HPC


This folder contains some scripts and templates to run OPFGenerator on a SLURM HPC cluster.

The file `submit_jobs.jl` is the main entry point and will automatically submit dependent jobs to the cluster. It is configured by default to use Pascal's AI4OPT charge account and the `embers` queue. You can change these settings by editing the first few lines of `submit_jobs.jl`.

To launch a dataset generation job, first create a TOML file with the desired parameters. Then, run the following command:

```bash
julia slurm/submit_jobs.jl <path/to/config/file.toml>
```