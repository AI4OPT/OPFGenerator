Copyright Georgia Tech 2022

[![Build][build-img]][build-url]

[build-img]: https://github.com/ai4opt/OPFGenerator/workflows/CI/badge.svg?branch=main
[build-url]: https://github.com/ai4opt/OPFGenerator/actions?query=workflow%3ACI

# OPFGenerator
Instance generator for various OPF problems (ACOPF & DCOPF currently supported)

- [OPFGenerator](#OPFGenerator)
  - [Installation instructions](#installation-instructions)
    - [Using HSL solvers](#using-hsl-solvers-ma27-ma57)
  - [Quick start](#quick-start)
  - [Generating datasets](#generating-datasets)
  - [Solution format](#solution-format)
    - [Primal variables](#primal-variables)
    - [Dual variables](#dual-variables)
  - [Datasets](#datasets)
    - [Format](#format)
    - [Loading from julia](#loading-from-julia)
    - [Loading from python](#loading-from-python)
  - [Loading and saving JSON files](#loading-and-saving-json-files)

## Installation instructions

This repository is a non-registered Julia package.

* Option 1: install as a Julia package. You can use it, but not modify the code
    ```julia
    using Pkg
    Pkg.add("git@github.com:AI4OPT/OPFGenerator.git")
    ```

* Option 2: clone the repository. Use this if you want to change the code
    ```bash
    git clone git@github.com:AI4OPT/OPFGenerator.git
    ```
    To use the package after cloning the repo
    ```bash
    $ cd OPFGenerator
    $ julia --project=.
    julia> using OPFGenerator
    ```

    If you are modifying the source code, it is recommened to use the package [`Revise.jl`](https://github.com/timholy/Revise.jl)
    so that you can use the changes without having to start Julia.
    Make sure you load `Revise` before loading `OPFGenerator` in your julia session.
    ```julia
    using Revise
    using OPFGenerator
    ```

### Using HSL solvers (ma27, ma57)

HSL solvers are supported via the [HSL.jl](https://github.com/JuliaSmoothOptimizers/HSL.jl) package (only v0.3.6 and above).
Follow installation instructions there.

To use HSL linear solvers with Ipopt, use the following:
```julia
using HSL
const LIB_COINHSL = HSL.libcoinhsl

using JuMP

solver = optimizer_with_attributes(Ipopt.Optimizer,
    "hsllib" => LIB_COINHSL,
    "linear_solver" => "ma27"  # or solver of your choice
)
acopf = OPFGenerator.build_acopf(data, solver)
```

## Quick start

```julia
using Random, PGLib, PowerModels
using OPFGenerator
PowerModels.silence()

using StableRNGs
rng = StableRNG(42)

old_data = make_basic_network(pglib("3_lmbd"))

# Load scaler using global scaling + uncorrelated LogNormal noise
config = Dict(
    "load" => Dict(
        "noise_type" => "ScaledLogNormal",
        "l" => 0.8,
        "u" => 1.2,
        "sigma" => 0.05,        
    )
)
opf_sampler  = SimpleOPFSampler(old_data, config)

# Generate a new instance
new_data = rand(rng, opf_sampler)

old_data["load"]["1"]["pd"]  # old 
1.1

new_data["load"]["1"]["pd"]  # new
1.1596456429775048
```

To generate multiple instances, run the above code in a loop
```julia
dataset = [
    OPFGenerator.rand(rng, opf_sampler)
    for i in 1:100
]
```

## Generating datasets

A script for generating multiple ACOPF instances is given in [`exp/sampler.jl`](exp/sampler.jl).

It is called from the command-line as follows:
```bash
julia --project=. exp/sampler.jl <path/to/config.toml> <seed_min> <seed_max>
```
where
* `<path/to/config.toml>` is a path to a valid configuration file in TOML format (see [`exp/ieee300.toml`](exp/ieee300.toml) for an example)
* `<seed_min>` and `<seed_max>` are minimum and maximum values for the random seed. Must be integer values.
    The script will generate instances for values `smin, smin+1, ..., smax-1, smax`.

## Solution format

The solution dictionary follows the [PowerModels result data format](https://lanl-ansi.github.io/PowerModels.jl/stable/result-data/).

### ACOPF Primal variables

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

### ACOPF Dual variables

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
|           | `"mu_va_diff"` | Voltage angle difference
| generator | `"mu_pg_lb"` | Active power generation lower bound
|           | `"mu_pg_ub"` | Active power generation upper bound
|           | `"mu_qg_lb"` | Reactive power generation lower bound
|           | `"mu_qg_ub"` | Reactive power generation upper bound

### DCOPF Primal variables

| Component | Key | Description |
|:---------:|:----|:------------|
| bus       | `"va"` | Voltage angle
| generator | `"pg"` | Power generation
| branch    | `"pf"` | Power flow

### DCOPF Dual variables

| Component | Key | Constraint  |
|:---------:|:----|:------------|
| bus       | `"lam_kirchhoff"` | Power balance
| generator | `"mu_pg_lb"` | Power generation lower bound
|           | `"mu_pg_ub"` | Power generation upper bound
| branch    | `"lam_ohm"` | Ohm's law
|           | `"mu_va_diff"` | Voltage angle difference

## Datasets

### Format

Each dataset is stored in an `.h5` file, organized as follows:
```
/
|-- meta
    |-- ref
    |-- seed
|-- input
    |-- pd
    |-- qd
    |-- br_status
|-- acopf_solution
    |-- meta
        |-- termination_status
        |-- primal_status
        |-- dual_status
        |-- solve_time
    |-- primal
        |-- vm
        |-- va
        |-- pg
        |-- qg
        |-- pf
        |-- qf
        |-- pt
        |-- qt
    |-- dual
        |-- mu_vm_lb
        |-- mu_vm_ub
        |-- lam_kircchoff_active
        |-- lam_kircchoff_reactive
        |-- mu_pg_lb
        |-- mu_pg_ub
        |-- mu_qg_lb
        |-- mu_qg_ub
        |-- mu_sm_fr
        |-- mu_sm_to
        |-- lam_ohm_active_fr
        |-- lam_ohm_active_to
        |-- lam_ohm_reactive_fr
        |-- lam_ohm_reactive_to
        |-- mu_va_diff
|-- dcopf_solution
    |-- meta
        |-- termination_status
        |-- primal_status
        |-- dual_status
        |-- solve_time
    |-- primal
        |-- va
        |-- pg
        |-- pf
    |-- dual
        |-- lam_kirchhoff
        |-- mu_pg_lb
        |-- mu_pg_ub
        |-- lam_ohm
        |-- mu_va_diff

```

### Loading from julia

```julia
using HDF5

D = h5read("dataset.h5", "/")  # read all the dataset into a dictionary
```

### Loading from python

The following code provides a starting point to load h5 datasets in python
```py

import h5py
import numpy as np


def parse_hdf5(path: str, preserve_shape: bool=False):
    dat = dict()

    def read_direct(dataset: h5py.Dataset):
        arr = np.empty(dataset.shape, dtype=dataset.dtype)
        dataset.read_direct(arr)
        return arr

    def store(name: str, obj):
        if isinstance(obj, h5py.Group):
            return
        elif isinstance(obj, h5py.Dataset):
            dat[name] = read_direct(obj)
        else:
            raise ValueError(f"Unexcepted type: {type(obj)} under name {name}. Expected h5py.Group or h5py.Dataset.")

    with h5py.File(path, "r") as f:
        f.visititems(store)

    if preserve_shape:
        # recursively correct the shape of the dictionary
        ret = dict()

        def r_correct_shape(d: dict, ret: dict):
            for k in list(d.keys()):
                if "/" in k:
                    k1, k2 = k.split("/", 1)
                    if k1 not in ret:
                        ret[k1] = dict()
                    r_correct_shape({k2: d[k]}, ret[k1])
                    del d[k]
                else:
                    ret[k] = d[k]

        r_correct_shape(dat, ret)

        return ret
    else:
        return dat
```

## Loading and saving JSON files

Use the `load_json` and `save_json` functions to load/save data to/from JSON files.
Uncompressed (`.json`) and compressed (`.json.gz` and `.json.bz2`) are supported automatically.

```julia
using OPFGenerator

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
