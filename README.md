Copyright Georgia Tech 2022-2024

[![Build][build-img]][build-url]
[![codecov][codecov-img]][codecov-url]
 
[build-img]: https://github.com/ai4opt/OPFGenerator/workflows/CI/badge.svg?branch=main
[build-url]: https://github.com/ai4opt/OPFGenerator/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/AI4OPT/OPFGenerator/graph/badge.svg
[codecov-url]: https://codecov.io/gh/AI4OPT/OPFGenerator

# OPFGenerator
Instance generator for various OPF problems (ACOPF & DCOPF currently supported)

- [OPFGenerator](#opfgenerator)
  - [Installation instructions](#installation-instructions)
    - [Using HSL solvers (ma27, ma57)](#using-hsl-solvers-ma27-ma57)
  - [Quick start](#quick-start)
    - [Generating random instances](#generating-random-instances)
    - [Building and solving OPF problems](#building-and-solving-opf-problems)
  - [Generating datasets](#generating-datasets)
  - [Solution format](#solution-format)
    - [ACOPF](#ACOPF)
    - [DCOPF](#DCOPF)
    - [SOCWRPowerModel](#socwrpowermodel)
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

Please refer to [`Ipopt.jl` installation instructions](https://github.com/jump-dev/Ipopt.jl?tab=readme-ov-file#linear-solvers)
    for how to install non-default linear solvers, e.g., HSL, Pardiso, etc...
Note that a specific license may be required.

To use HSL linear solvers when solving OPF instances, set the parameter "linear_solver" to "ma27" or "ma57" in the config file.
The recommended solver for Ipopt is `ma27`.
```toml
solver.name = "Ipopt"
solver.attributes.linear_solver = "ma27"
```

## Quick start

### Building and solving an OPF problem

```julia
using PGLib, PowerModels
using OPFGenerator
using Ipopt

# Load data
pm_network = PowerModels.make_basic_network(pglib("14_ieee"))
data = OPFGenerator.OPFData(pm_network)

# Build OPF optimization model
# To switch formulations, replace ACOPF with, e.g., DCOPF or SOCOPF
opf = OPFGenerator.build_opf(OPFGenerator.ACOPF, data, Ipopt.Optimizer)

# Solve and extract result
OPFGenerator.solve!(opf)
res = OPFGenerator.extract_result(opf)
```



### Generating random OPF instances

```julia
using Random, PGLib, PowerModels
using OPFGenerator
PowerModels.silence()

rng = MersenneTwister(42)

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
* `<path/to/config.toml>` is a path to a valid configuration file in TOML format (see [`exp/config.toml`](exp/config.toml) for an example)
* `<seed_min>` and `<seed_max>` are minimum and maximum values for the random seed. Must be integer values.
    The script will generate instances for values `smin, smin+1, ..., smax-1, smax`.

## Solution format


Individual solutions are stored in a dictionary that follows the [PowerModels result data format](https://lanl-ansi.github.io/PowerModels.jl/stable/result-data/).
As a convention, dual variables of equality constraints are named `lam_<constraint_ref>`, dual variables of (scalar) inequality constraints are named `mu_<constraint_ref>`, and dual variables of conic constraints are named `nu_<constraint_ref>_<i>` where `i` denotes the coordinate index.

Collections of instances & solutions are stored in a compact, array-based HDF5 file (see [Datasets/Format](#format) below.)

### ACOPF

See [PowerModels documentation](https://lanl-ansi.github.io/PowerModels.jl/stable/formulation-details/#PowerModels.ACPPowerModel).

Primal variables

| Component | Key | Description |
|:---------:|:----|:------------|
| bus       | `"vm"` | Nodal voltage magnitude
|           | `"va"` | Nodal voltage angle
| generator | `"pg"` | Active power generation
|           | `"qg"` | Reactive power generation
| branch    | `"pf"` | Branch active power flow (fr)
|           | `"pt"` | Branch active power flow (to)
|           | `"qf"` | Branch reactive power flow (fr)
|           | `"qt"` | Branch reactive power flow (to)

Dual variables

| Component | Key                        | Constraint |
|:---------:|:---------------------------|:------------|
| bus       | `"mu_vm_lb"`               | Nodal voltage magnitude lower bound
|           | `"mu_vm_ub"`               | Nodal voltage magnitude upper bound
|           | `"lam_kirchhoff_active"`   | Nodal active power balance
|           | `"lam_kirchhoff_reactive"` | Nodal reactive power balance
| generator | `"mu_pg_lb"`               | Active power generation lower bound
|           | `"mu_pg_ub"`               | Active power generation upper bound
|           | `"mu_qg_lb"`               | Reactive power generation lower bound
|           | `"mu_qg_ub"`               | Reactive power generation upper bound
| branch    | `"mu_sm_fr"`               | Thermal limit (fr)
|           | `"mu_sm_to"`               | Thermal limit (to)
|           | `"lam_ohm_active_fr"`      | Ohm's law; active power (fr)
|           | `"lam_ohm_active_to"`      | Ohm's law; active power (to)
|           | `"lam_ohm_reactive_fr"`    | Ohm's law; reactive power (fr)
|           | `"lam_ohm_reactive_to"`    | Ohm's law; reactive power (to)
|           | `"mu_va_diff"`             | Voltage angle difference

### DCOPF

See [PowerModels documentation](https://lanl-ansi.github.io/PowerModels.jl/stable/formulation-details/#PowerModels.DCPPowerModel).

Primal variables

| Component | Key | Description |
|:---------:|:----|:------------|
| bus       | `"va"` | Voltage angle
| generator | `"pg"` | Power generation
| branch    | `"pf"` | Power flow

Dual variables

| Component | Key               | Constraint  |
|:---------:|:------------------|:------------|
| bus       | `"lam_kirchhoff"` | Power balance
| generator | `"mu_pg_lb"`      | Power generation lower bound
|           | `"mu_pg_ub"`      | Power generation upper bound
| branch    | `"lam_ohm"`       | Ohm's law
|           | `"mu_sm_lb"`      | Thermal limit (lower bound)
|           | `"mu_sm_ub"`      | Thermal limit (upper bound)
|           | `"mu_va_diff"`    | Voltage angle difference

### SOCWRPowerModel

See [PowerModels documentation](https://lanl-ansi.github.io/PowerModels.jl/stable/formulation-details/#PowerModels.SOCWRPowerModel).

Primal variables

| Component | Key | Description |
|:---------:|:----|:------------|
| bus       | `"w"`  | Squared voltage magnitude
| generator | `"pg"` | Active power generation
|           | `"qg"` | Reactive power generation
| branch    | `"pf"` | Branch active power flow (fr)
|           | `"pt"` | Branch active power flow (to)
|           | `"qf"` | Branch reactive power flow (fr)
|           | `"qt"` | Branch reactive power flow (to)
|           | `"wr"` | Real part of voltage product
|           | `"wi"` | Imaginary part of voltage product

Dual variables depend on whether a quadratic (`SOCOPFQuad`) or conic (`SOCOPF`) formulation is considered.
The former are identified with `(quad)`, the latter with `(cone)` in the table below.
Only the dual variables of quadratic constraints are affected by this distinction.
Note that conic dual are high-dimensional variables, and separate coordinates are stored separately.

| Component | Key                        | Constraint |
|:---------:|:---------------------------|:------------|
| bus       | `"mu_vm_lb"`               | Nodal voltage magnitude lower bound
|           | `"mu_vm_ub"`               | Nodal voltage magnitude upper bound
|           | `"lam_kirchhoff_active"`   | Nodal active power balance
|           | `"lam_kirchhoff_reactive"` | Nodal reactive power balance
| generator | `"mu_pg_lb"`               | Active power generation lower bound
|           | `"mu_pg_ub"`               | Active power generation upper bound
|           | `"mu_qg_lb"`               | Reactive power generation lower bound
|           | `"mu_qg_ub"`               | Reactive power generation upper bound
| branch    | `"lam_ohm_active_fr"`      | Ohm's law; active power (fr)
|           | `"lam_ohm_active_to"`      | Ohm's law; active power (to)
|           | `"lam_ohm_reactive_fr"`    | Ohm's law; reactive power (fr)
|           | `"lam_ohm_reactive_to"`    | Ohm's law; reactive power (to)
|           | `"mu_va_diff_lb"`          | Voltage angle difference lower bound
|           | `"mu_va_diff_ub"`          | Voltage angle difference upper bound
|           | `"mu_sm_fr"`               | (quad) Thermal limit (fr)
|           | `"mu_sm_to"`               | (quad) Thermal limit (to)
|           | `"mu_voltage_prod_quad"`   | (quad) Voltage product relaxation (Jabr)
|           | `"nu_voltage_prod_soc_1"`  | (cone) Voltage product relaxation (Jabr)
|           | `"nu_voltage_prod_soc_2"`  | (cone) Voltage product relaxation (Jabr)
|           | `"nu_voltage_prod_soc_3"`  | (cone) Voltage product relaxation (Jabr)
|           | `"nu_voltage_prod_soc_4"`  | (cone) Voltage product relaxation (Jabr)
|           | `"nu_sm_fr_1"`             | (cone) Thermal limit (fr)
|           | `"nu_sm_fr_2"`             | (cone) Thermal limit (fr)
|           | `"nu_sm_fr_3"`             | (cone) Thermal limit (fr)
|           | `"nu_sm_to_1"`             | (cone) Thermal limit (to)
|           | `"nu_sm_to_2"`             | (cone) Thermal limit (to)
|           | `"nu_sm_to_3"`             | (cone) Thermal limit (to)


## Datasets

### Format

Each dataset is stored in an `.h5` file, organized as follows.
See [Solution format](#solution-format) for a list of each formulation's primal and dual variables.
```
/
|-- meta
|   |-- ref
|   |-- config
|-- input
|   |-- seed
|   |-- pd
|   |-- qd
|   |-- br_status
|   |-- gen_status
|-- ACOPF
|   |-- meta
|   |   |-- termination_status
|   |   |-- primal_status
|   |   |-- dual_status
|   |   |-- solve_time
|   |-- primal
|   |   |-- pg
|   |   |-- qg
|   |   |-- ...
|   |-- dual
|       |-- lam_kirchhoff_active
|       |-- lam_kirchhoff_reactive
|       |-- ...
|-- DCOPF
|   |-- meta
|   |-- primal
|   |-- dual
|-- SOCOPF
    |-- meta
    |-- primal
    |-- dual
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

## Acknowledgements

This material is based upon work supported by the National Science Foundation under Grant No. 2112533 NSF AI Institute for Advances in Optimization ([AI4OPT](https://www.ai4opt.org/)). 
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
