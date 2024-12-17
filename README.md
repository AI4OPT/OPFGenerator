Copyright Georgia Tech 2022-2024

[![Build][build-img]][build-url]
[![codecov][codecov-img]][codecov-url]
 
[build-img]: https://github.com/ai4opt/OPFGenerator/workflows/CI/badge.svg?branch=main
[build-url]: https://github.com/ai4opt/OPFGenerator/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/AI4OPT/OPFGenerator/graph/badge.svg
[codecov-url]: https://codecov.io/gh/AI4OPT/OPFGenerator

# OPFGenerator
Instance generator for various OPF problems (ACOPF & DCOPF currently supported)

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

# Load data from PGLib case
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
using Random 
using PGLib, PowerModels
using OPFGenerator

# Load data from PGLib case
pm_network = PowerModels.make_basic_network(pglib("14_ieee"))
data = OPFGenerator.OPFData(pm_network)

# Load scaler using global scaling + uncorrelated LogNormal noise
config = Dict(
    "load" => Dict(
        "noise_type" => "ScaledUniform",
        "l" => 0.8,
        "u" => 1.2,
        "sigma" => 0.10,       
    )
)
opf_sampler = SimpleOPFSampler(data, config)

# Generate a new instance
rng = MersenneTwister(42)
new_data = rand(rng, opf_sampler)

data.pd[1] # old 
0.217

new_data.pd[1]  # new
0.21480423013573954  # ⚠️exact value depends on Julia version)
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

## Datasets

### File structure

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

See each formulation's documentation for details on the format of each `.h5` file.

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
