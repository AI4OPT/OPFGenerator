"""
    convert_float_data(D, F)

Convert all floating-point scalars and arrays to `F`.

# Arguments
* `D`: Should be a JSON-serializable dictionary, which generally means that
    all keys are `String` and all values are JSON-compatible.
* `F`: Must be a subtype of `AbstractFloat`

# Returns
* `d::Dict{String,Any}`: a dictionary with same nested structure as `D`,
    with all floating-point scalars and arrays converted to `F`.
"""
function convert_float_data(D, F)
    F <: AbstractFloat || error("`F` must be a subtype of `AbstractFloat`")
    keytype(D) <: AbstractString || error("Only `AbstractString`-like keys are supported.")
    
    d = Dict{String,Any}()

    for (k, v) in D
        if isa(v, Dict)
            d[k] = convert_float_data(v, F)
        elseif isa(v, AbstractFloat)
            d[k] = convert(F, v)
        elseif isa(v, Array{T,N} where {T <: AbstractFloat, N})
            d[k] = convert(Array{F,ndims(v)}, v)
        else
            d[k] = v
        end
    end

    return d
end
