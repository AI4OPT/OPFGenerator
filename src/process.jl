using Base.Iterators
using Base.Threads
using ProgressMeter
using HDF5

"""
    initialize_res(data)

Create an empty dataset for original instance `data`.

This function outputs a dictionary `D` that is meant to be passed to
    [`add_datapoint!`](@ref).
"""
function initialize_res(config)
    casename = config["ref"]
    data     = make_basic_network(pglib(casename))

    # Grab some information
    opf_formulations = sort(collect(keys(config["OPF"])))

    # initialize dataset
    D = Dict{String,Any}(
        "meta" => Dict{String,Any}(
            "ref"  => data["name"],
            "config" => deepcopy(config),
        ),
        "input" => Dict{String,Any}(
            "seed"      => Int[],
            "pd"        => Vector{Float64}[],
            "qd"        => Vector{Float64}[],
            "br_status" => Vector{Bool}[],
        ),
    )

    for opf_formulation in opf_formulations
        D[opf_formulation] = Dict{String,Any}(
            "meta" => Dict{String,Any}(
                "termination_status" => String[],
                "primal_status" => String[],
                "dual_status" => String[],
                "solve_time" => Float64[],
            ),
            "primal" => Dict{String,Any}(),
            "dual"   => Dict{String,Any}(),
        )
    end

    return D
end

function add_datapoint!(D, d)
    config = D["meta"]["config"]

    D_input::Dict{String,Any}    = D["input"]
    meta::Dict{String,Any} = d["meta"]
    data::Dict{String,Any} = d["data"]

    N = length(data["bus"])
    E = length(data["branch"])
    L = length(data["load"])
    G = length(data["gen"])

    # TODO: sanity checks
    @assert d["meta"]["ref"] == config["ref"]
    @assert all(haskey(d, k) for (k, v) in config["OPF"])

    # Meta info
    seed::Int = meta["seed"]
    push!(D["input"]["seed"], seed)

    # Input data
    pd::Vector{Float64} = [data["load"]["$l"]["pd"] for l in 1:L]
    qd::Vector{Float64} = [data["load"]["$l"]["qd"] for l in 1:L]
    br::Vector{Bool}    = [Bool(data["branch"]["$e"]["br_status"]) for e in 1:E]
    push!(D_input["pd"], pd)
    push!(D_input["qd"], qd)
    push!(D_input["br_status"], br)

    # Add OPF solutions
    opf_formulations = sort(collect(keys(config["OPF"])))
    for opf_formulation in opf_formulations
        Dsol = get!(D, opf_formulation, Dict{String,Any}())
        resh5 = json2h5(OPF2TYPE[opf_formulation], d[opf_formulation])
        for cat in keys(resh5)
            get!(Dsol, cat, Dict{String,Any}())
            for (k, v) in resh5[cat]
                if !haskey(Dsol[cat], k)
                    Dsol[cat][k] = [v]  # Make sure that all relevant keys exist
                else
                    push!(Dsol[cat][k], v)
                end
            end
        end
    end

    return D
end

function save_h5(filename::AbstractString, D)
    h5open(filename, "w") do file
        save_h5(file, D)
    end
    return nothing
end

function save_h5(file::HDF5.File, D::Dict)
    for (k, v) in D
        @assert isa(v, Dict)
        gr = create_group(file, k)
        save_h5(gr, v)
    end
    return nothing
end

function save_h5(gr::HDF5.Group, D::Dict)
    for (k, v) in D
        if isa(v, Array)
            gr[k] = v
        elseif isa(v, AbstractString)
            gr[k] = string(v)
        elseif isa(v, Union{Int64,Float64})
            gr[k] = v
        elseif isa(v, Dict)
            gr_ = create_group(gr, k)
            save_h5(gr_, v)
        else
            error("Unsupported data type for writing to a group: $k::$(typeof(v))")
        end
    end
    return nothing
end

"""
    convert_to_h5!(D)

Convert dataset `D` from list-of-list to h5-compatible format.

The input dictionary `D` should be the output of [`initialize_res`](@ref),
    to which an arbitrary number of datapoints have been added.
This function replaces non h5-compatible arrays (namely vectors of vectors)
    with contiguous arrays (`Matrix`).
The dataset `D` is modified in-place.
"""
function convert_to_h5!(D::Dict)
    config = D["meta"]["config"]
    # Input data
    dat  = D["input"]
    for k in ["pd", "qd", "br_status"]
        dat[k] = _vecvec2mat(dat[k])
    end
        
    opf_formulations = sort(collect(keys(config["OPF"])))
    for opf_formulation in opf_formulations
        opfres = D[opf_formulation]
        for cat in ["meta", "primal", "dual"]
            for (k, v) in opfres[cat]
                opfres[cat][k] = _vecvec2mat(v)
            end
        end
    end

    return nothing
end

_vecvec2mat(V) = reduce(hcat, V)

function parse_jsons(config::Dict;
    show_progress::Bool=true,
    batch_size::Int=0,
    force_process::Bool=true,
)
    exp_folder = config["export_dir"]
    json_folder = joinpath(exp_folder, "res_json")
    h5_folder   = joinpath(exp_folder, "res_h5")

    all_files = filter(s -> endswith(s, r".json|.json.gz"), readdir(json_folder))
    sort!(all_files)
    N = length(all_files)

    # Read first file
    fname = all_files[1]
    ref = load_json(joinpath(json_folder, fname))["meta"]["ref"]
    data = make_basic_network(pglib(ref))

    # If no batch size was provided, process everything in one go
    batch_size = (batch_size <= 0) ? N : batch_size
    F = partition(all_files, batch_size)
    println("Processing $N results files from $(json_folder)")
    println("Chunk size: $(batch_size) ($(length(F)) chunks)")
    println("Results will be exported to $(joinpath(h5_folder, ref * "<chunk_index>.h5."))")

    p = Progress(N; enabled=show_progress)

    for (b, files) in enumerate(F)
        f5name = joinpath(h5_folder, ref * "_$(b).h5")
        !force_process && isfile(f5name) && continue

        D = initialize_res(config)

        # Process all files by chunks
        L = ReentrantLock()
        num_read_error = Threads.Atomic{Int}(0)
        @threads for fname in files
            next!(p)
            local d = try
                load_json(joinpath(json_folder, fname))
            catch err
                @info "Error while reading $(fname)" err
                num_read_error[] += 1
                continue
            end
            lock(L) do
                add_datapoint!(D, d)
            end
        end

        (num_read_error[] > 0) && @info "$(num_read_error[]) errors importing JSON files while processing batch $b."
        @info "Saving batch $b to disk"
        convert_to_h5!(D)
        save_h5(f5name, D)
        GC.gc()  # helps keep memory under control
    end

    return nothing
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

function _merge_h5!(D, args...)
    N = length(args)
    all(d -> isa(d, AbstractDict), args) || throw(ArgumentError("All arguments must be dictionaries"))
    for (k, v) in D
        if isa(v, AbstractVector)
            # append both vectors
            D[k] = reduce(vcat, [d[k] for d in args])
        elseif isa(v, AbstractMatrix)
            # concatenate both matrices
            D[k] = reduce(hcat, [d[k] for d in args])
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

The dictionary `D` should be in h5-compatible format.
It is modified in-place.
"""
function _sort_h5!(D::Dict{String,Any})
    p = sortperm(D["input"]["seed"])

    _sort_h5!(D, p)
    return nothing
end

function _sort_h5!(D::Dict{String,Any}, p::Vector{Int})
    for (k, v) in D
        if isa(v, AbstractVector)
            D[k] = v[p]
        elseif isa(v, AbstractMatrix)
            D[k] = v[:, p]
        elseif isa(v, AbstractDict)
            _sort_h5!(v, p)
        end
    end
    return D
end
