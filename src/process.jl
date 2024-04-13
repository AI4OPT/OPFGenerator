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
function tensorize(V::Vector{T}) where {T <: Union{String,Number}}
    length(V) > 0 || error("Trying to tensorize an empty collection")
    return V
end

function tensorize(V::Vector{Array{T,N}}) where {T,N}
    length(V) > 0 || error("Trying to tensorize an empty collection")
    return stack(V)
end

function _merge_h5(V::Vector{<:Dict})
    length(V) > 0 || error("Cannot merge an empty collection")

    v0 = V[1]
    D = Dict{String,Any}()
    for k in keys(v0)
        D[k] = _merge_h5([v[k] for v in V])
    end

    return D
end

function _merge_h5(V::Vector{Array{T,N}}) where{T,N}
    return cat(V..., dims=ndims(first(V)))
end

"""
    _sort_h5!(D)

Sort dataset `D` in increasing order of random seeds.

The dictionary `D` should be in h5-compatible format. It is modified in-place.

The function expects `D["seed"]` or `D["meta"]["seed"]` to exist and be an array that can sorted.
    An error is thrown if this entry is not found.
"""
function _sort_h5!(D::Dict{String,Any})
    haskey(D, "meta") && haskey(D["meta"], "seed") || error("Invalid H5 dataset: missing random seeds in the \"meta\" section.")
    seeds::Vector{Int} = D["meta"]["seed"]
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

function _dedupe!(D::Dict{String,Any})
    seeds = D["meta"]["seed"]
    unique_seeds_idx = unique(i -> seeds[i], eachindex(seeds))
    if length(unique_seeds_idx) == length(seeds)
        return D
    end

    _dedupe!(D, unique_seeds_idx)
end

function _dedupe!(D::Dict{String,Any}, unique_seeds_idx::Vector{Int})
    for (k, v) in D
        if isa(v, Array)
            D[k] = convert(typeof(v), selectdim(v, ndims(v), unique_seeds_idx))
        elseif isa(v, AbstractDict)
            _dedupe!(v, unique_seeds_idx)
        elseif isa(v, Union{String,Number})
            # nothing to do
        else
            error("Unexpected type $(typeof(v)) for entry $k while deduping H5 dataset")
        end
    end
    return D
end