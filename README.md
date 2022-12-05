Copyright Georgia Tech 2022

# ACOPFGenerator
Instance generator for ACOPF problem

- [ACOPFGenerator](#acopfgenerator)
  - [Installation instructions](#installation-instructions)
  - [Quick start](#quick-start)
  - [Solution format](#solution-format)
    - [Primal variables](#primal-variables)
    - [Dual variables](#dual-variables)
  - [Loading and saving JSON files](#loading-and-saving-json-files)

## Installation instructions

This repository is a non-registered Julia package.

* Option 1: install as a Julia package. You can use it, but not modify the code
    ```julia
    using Pkg
    Pkg.add("git@github.com:AI4OPT/ACOPFGenerator.git")
    ```

* Option 2: clone the repository. Use this if you want to change the code
    ```bash
    git clone git@github.com:AI4OPT/ACOPFGenerator.git
    ```
    To use the package after cloning the repo
    ```bash
    $ cd ACOPFGenerator
    $ julia --project=.
    julia> using ACOPFGenerator
    ```

    If you are modifying the source code, it is recommened to use the package [`Revise.jl`](https://github.com/timholy/Revise.jl)
    so that you can use the changes without having to start Julia.
    Make sure you load `Revise` before loading `ACOPFGenerator` in your julia session.
    ```julia
    using Revise
    using ACOPFGenerator
    ```

## Quick start

```julia
using Random, PGLib, PowerModels
using ACOPFGenerator
PowerModels.silence()

using StableRNGs
rng = StableRNG(42)

data = make_basic_network(pglib("3_lmbd"))

# Uncorrelated, uniformly-distributed multiplicative noise
ls_uniform = SimpleLoadScaling(data, 0.8, 1.2)
# Correlated scaling + uncorrelated noise
ls_lognorm = ScaleLogNorm(data, 0.8, 1.2, 0.01)

# Each sampler may be used 
opf_sampler_uniform  = SimpleOPFSampler(data, ls_uniform)
opf_sampler_lognorm  = SimpleOPFSampler(data, ls_lognorm)

# Generate a new instance. First we 
new_data_uniform = ACOPFGenerator.sample(rng, opf_sampler_uniform)
new_data_lognorm = ACOPFGenerator.sample(rng, opf_sampler_lognorm)

data["load"]["1"]["pd"]
1.1

new_data_uniform["load"]["1"]["pd"]
1.135426539581486

new_data_lognorm["load"]["1"]["pd"]
1.208580250500669
```

## Solution format

The solution dictionary follows the [PowerModels result data format](https://lanl-ansi.github.io/PowerModels.jl/stable/result-data/).

### Primal variables

| Component | Key | Description |
|:---------:|:----|:------------|
| bus       | `"vm"` | Nodal voltage magnitude
|           | `"va"` | Nodal voltage angle
| branch    | `"pf"` | Branch active power flow (fr)
|           | `"pt"` | Branch active power flow (to)
|           | `"qf"` | Branch reactive power flow (fr)
|           | `"qt"` | Branch reactive power flow (to)
| generator | `"pg"` | Active power generation
|           | `"qg"` | Reactive power generation

### Dual variables

As a convention, dual variables of equality constraints are named `lam_<constraint_ref>`, and dual variables of inequality constraints are named `mu_<constraint_ref>`.
Dual variables are stored together with each component's primal solution.

| Component | Key | Constraint |
|:---------:|:----|:------------|
| bus       | `"mu_vm_lb"` | Nodal voltage magnitude lower bound
|           | `"mu_vm_ub"` | Nodal voltage magnitude upper bound
|           | `"lam_pb_active"` | Nodal active power balance
|           | `"lam_pb_reactive"` | Nodal reactive power blaance
| branch    | `"mu_sm_fr"` | Thermal limit (fr)
|           | `"mu_sm_to"` | Thermal limit (to)
|           | `"lam_ohm_active_fr"` | Ohm's law; active power (fr)
|           | `"lam_ohm_active_to"` | Ohm's law; active power (to)
|           | `"lam_ohm_reactive_fr"` | Ohm's law; reactive power (fr)
|           | `"lam_ohm_reactive_to"` | Ohm's law; reactive power (to)
| generator | `"mu_pg_lb"` | Active power generation lower bound
|           | `"mu_pg_ub"` | Active power generation upper bound
|           | `"mu_qg_lb"` | Reactive power generation lower bound
|           | `"mu_qg_ub"` | Reactive power generation upper bound


## Loading and saving JSON files

Use the `load_json` and `save_json` functions to load/save data to/from JSON files.
Uncompressed (`.json`) and compressed (`.json.gz` and `.json.bz2`) are supported automatically.

```julia
using ACOPFGenerator

# Load a dictionary from a JSON file
d = load_json("my_json_file.json")
d = load_json("my_json_file.json.gz")
d = load_json("my_json_file.json.bz2")

# Save a dictionary to JSON file
save_json("my_new_jsonfile.json", d)
save_json("my_pretty_jsonfile.json", d, indent=2)  # prettier formatting
save_json("my_new_jsonfile.json.gz", d)
save_json("my_new_jsonfile.json.bz2", d)
```
