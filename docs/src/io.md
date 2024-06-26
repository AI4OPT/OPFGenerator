`OPFGenerator` provides utilities to save/load files in the HDF5 and JSON formats.

## HDF5

To load HDF5 files, use the [`load_h5`](@ref), which is the same as `HDF5.h5read`.
```julia
h = load_h5("test_file.h5", "/")  # load all HDF5 file
```

To save HDF5 files to disk, use the [`save_h5`](@ref) function.
All keys of the dictionary (and sub-dictionaries) must be `String`s.

```julia
save_h5("myfile.h5", d)
```

!!! warning
    
    When exporting data to an HDF5 file, only the following data types are supported:
    * `String`
    * Number types [supported by HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/#Supported-data-types)
    * `Array`s of those

    Numerical data in an unsupported type will be converted to `Float64`, when possible.

## JSON

!!! warning
    Only `.json` extensions are supported

Use the [`load_json`](@ref) and [`save_json`](@ref) functions to load/save data to/from JSON files.

```julia
using OPFGenerator

# Load a dictionary from a JSON file
d = load_json("my_json_file.json")

# Save a dictionary to a JSON file
save_json("my_new_jsonfile.json", d)
save_json("my_pretty_jsonfile.json", d, indent=2)  # prettier formatting
```

Compressed JSON files (`.json.gz` and `.json.bz2`) are supported automatically

```julia
# Load a dictionary from a compressed JSON file
d = load_json("my_json_file.json.gz")
d = load_json("my_json_file.json.bz2")

# Save a dictionary to a compressed JSON file
save_json("my_new_jsonfile.json.gz", d)
save_json("my_new_jsonfile.json.bz2", d)
```
