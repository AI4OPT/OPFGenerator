# ACOPFGenerator
Instance generator for ACOPF problem

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

