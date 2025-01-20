using HDF5

# Supported datatypes as per https://juliaio.github.io/HDF5.jl/v0.17/#Supported-data-types
const HDF5_SUPPORTED_NUMBER_TYPES = Union{
    Bool,
    UInt8, UInt16, UInt32, UInt64,
     Int8,  Int16,  Int32,  Int64,
    Float32, Float64,
    # Complex versions of all of these
    Complex{Bool},
    Complex{UInt8}, Complex{UInt16}, Complex{UInt32}, Complex{UInt64},
    Complex{Int8},  Complex{Int16},  Complex{Int32},  Complex{Int64},
    Complex{Float32}, Complex{Float64},
}

"""
    load_h5
"""
load_h5 = HDF5.h5read

"""
    save_h5(filename, D; warn=true)

Saves dictionary `D` to HDF5 file `filename`.

# Arguments
* `filename::AbstractString`: Path to the HDF5 file; must be a valid path.
* `D`: Dictionary to save to the file. 
    All keys in `D` must be of `String` type, and it must be HDF5-compatible.
    Additional restrictions are enforced on the values of `D`, see below.
* `warn::Bool=true`: Whether to raise a warning when converting numerical data.

!!! warning
    Only the following types are supported:
    * `String`
    * (un)signed integers up to 64-bit precision
    * `Float32` and `Float64`
    * `Complex` versions of the above numeric types
    * Dense `Array`s of the the above scalar types

    Numerical data whose type is not listed above will be converted to `Float64`.
    Set `warn=false` to disable such warnings.
    An error will be thrown if this conversion is not possible 
    
"""
function save_h5(filename::AbstractString, D; warn=true)
    h5open(filename, "w") do file
        _save_h5(file, D; warn)
    end
    return nothing
end

# Conversion for concrete types
_convert_to_h5_supported(::Type{T}) where{T<:AbstractString} = String
_convert_to_h5_supported(::Type{T}) where{T<:HDF5_SUPPORTED_NUMBER_TYPES} = T
_convert_to_h5_supported(::Type{T}) where{T<:Real} = Float64
_convert_to_h5_supported(::Type{T}) where{T<:Complex} = Complex{Float64}
_convert_to_h5_supported(::Type{T}) where {T} = Any
_convert_to_h5_supported(::Type{Array{T,N}}) where {T,N} = Array{_convert_to_h5_supported(T),N}

function _save_h5(fid::HDF5.H5DataStore, D::Dict; warn=true)
    for (k, v) in D
        if isa(v, Dict)
            # TODO: check that `k` is a `String` and does not contain any `/` characters
            # (these will create subgroups in the HDF5 file)
            gr = create_group(fid, k)
            _save_h5(gr, v)
        else
            T = typeof(v)
            T5 = _convert_to_h5_supported(T)
            # We could not convert to a supported type --> throw an error
            ((T5 == Any) || ((T5 <: AbstractArray) && eltype(T) == Any)) && error(
                """Unsupported data type for writing to an HDF5 group: \"$k\"::$(typeof(v)).
                HDF5 only supports the following types:
                * `String`
                * (un)signed integers up to 64-bit precision
                * `Float32` and `Float64`
                * `Complex` versions of the above numeric types
                * Dense `Array`s of the the above scalar types

                (see https://juliaio.github.io/HDF5.jl/v0.17/#Supported-data-types)

                Consider either converting to a different type, or using a different file format.\n"""
            )
            # Raise a warning if we are converting from a non-HDF5-supported numerical type
            warn && (T5 != T) && !(T <: AbstractString) && @warn(
                """Unsupported data type for writing to an HDF5 group: \"$k\"::$(T).
                This value was converted to $(T5), which may incur a loss of precision.""",
                maxlog=10,
            )
            fid[k] = convert(T5, v)
        end
    end
    return nothing
end
