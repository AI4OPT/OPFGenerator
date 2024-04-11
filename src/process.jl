using Base.Iterators
using Base.Threads

"""
    tensorize(V)

Concatenate elements of `V` into a higher-dimensional tensor.

Similar to `Base.stack`, with one major difference: if `V` is a vector of scalars,
    the result is a 2D array `M` whose last dimension is `length(V)`,
    and such that `M[:, i] == V[i]`.

This function is only defined for `Vector{T}` and `Vector{Array{T,N}}` inputs,
    to avoid any unexpected behavior of `Base.stack`.
"""
function tensorize(V::Vector{T}) where {T}
    # Check that all elements have same size
    length(V) > 0 || error("Trying to tensorize an empty collection")

    # We do not use `reduce(hcat, V)` to avoid type instability.
    # Since `V` may be an arbitrary collection, we explictly allocate the output
    return reshape(copy(V), (1, length(V)))
end

function tensorize(V::Vector{Array{T,N}}) where {T,N}
    length(V) > 0 || error("Trying to tensorize an empty collection")
    return stack(V)
end

function _merge_h5(args...)
    N = length(args)

    # Check that all arguments are Dict
    all(d -> isa(d, AbstractDict), args) || throw(ArgumentError("All arguments must be dictionaries"))

    N == 0 && return Dict{String,Any}()
    D = deepcopy(first(args))
    _merge_h5!(D, args...)

    return D
end

function _merge_h5_new(V::Vector{<:Dict})
    length(V) > 0 || error("Cannot merge an empty collection")

    v0 = V[1]
    D = Dict{String,Any}()
    for k in keys(v0)
        D[k] = _merge_h5_new([v[k] for v in V])
    end

    return D
end

function _merge_h5_new(V::Vector{Array{T,N}}) where{T,N}
    M = stack(V)
    return reshape(M, (size(M)[1:end-2]..., prod(size(M)[end-1:end])))
end

function _merge_h5!(D, args...)
    N = length(args)
    all(d -> isa(d, AbstractDict), args) || throw(ArgumentError("All arguments must be dictionaries"))
    for (k, v) in D
        if isa(v, Array)
            N = ndims(v)
            M = stack(d[k] for d in args)
            D[k] = reshape(M, (size(M)[1:end-2]..., prod(size(M)[end-1:end])))
        elseif isa(v, AbstractDict)
            # recursively merge sub-dictionaries
            _merge_h5!(D[k], [d[k] for d in args]...)
        else
            # Check that all values are the same
            all(d[k] == v for d in args) || error("Different values for entry $k of type $(typeof(v))")
        end
    end

    return nothing
end

"""
    _sort_h5!(D)

Sort dataset `D` in increasing order of random seeds.

The dictionary `D` should be in h5-compatible format. It is modified in-place.

The function expects `D["seed"]` or `D["meta"]["seed"]` to exist and be an array that can sorted.
    An error is thrown if this entry is not found.
"""
function _sort_h5!(D::Dict{String,Any})
    seeds = if haskey(D, "seed") 
        D["seed"][:]
    elseif haskey(D, "meta") && haskey(D["meta"], "seed")
        D["meta"]["seed"][:]
    else
        # safeguards
        error("Cannot find random seeds")
    end

    p = sortperm(seeds)

    _sort_h5!(D, p)
    return nothing
end

function _sort_h5!(D::Dict{String,Any}, p::Vector{Int}; checkperm=true)
    (!checkperm || isperm(p)) || error("Input permutation is not a valid permutation.")
    for (k, v) in D
        if isa(v, Array)
            # Sanity checks
            ns = size(v)
            length(p) == ns[end] || error("Invalid permutation size: entry \"$(k)\" has $(ns[end]) entries but permutation has size $(length(p)).")

            # flatten all dimensions except the last one
            M = reshape(v, (prod(ns[1:end-1]), ns[end]))
            # sort (flattened) columns
            M = M[:, p]
            # reshape back to original dimensions
            D[k] = reshape(M, ns)
        elseif isa(v, AbstractDict)
            # Sort sub-dictionaries
            # We do not need to re-check that `p` is a valid permutation, since either
            #   1) `checkperm` was `false`, or `checkperm` was `true` and we already checked
            _sort_h5!(v, p; checkperm=false)
        elseif isa(v, Union{String,Number})
            # nothing to do
        else
            # safeguard for unexpected types
            error("Unexpected type $(typeof(v)) for entry $k while sorting H5 dataset")
        end
    end
    return D
end
