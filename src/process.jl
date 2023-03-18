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
function initialize_res(data)
    # initialize dataset
    D = Dict{String,Any}(
        "meta" => Dict{String,Any}(
            "seed" => Int[],
            "ref"  => data["name"],
        ),
        "input" => Dict{String,Any}(
            "pd"        => Vector{Float64}[],
            "qd"        => Vector{Float64}[],
            "br_status" => Vector{Bool}[],
        ),
        "solution" => Dict{String,Any}(
            "meta" => Dict{String,Any}(
                "termination_status" => String[],
                "primal_status" => String[],
                "dual_status" => String[],
                "solve_time" => Float64[],
            ),
            "primal" => Dict{String,Any}(
                "vm" => Vector{Float64}[],
                "va" => Vector{Float64}[],
                "pg" => Vector{Float64}[],
                "qg" => Vector{Float64}[],
                "pf" => Vector{Float64}[],
                "qf" => Vector{Float64}[],
                "pt" => Vector{Float64}[],
                "qt" => Vector{Float64}[],
            ),
            "dual" => Dict{String,Any}(
                "mu_vm_lb" => Vector{Float64}[],
                "mu_vm_ub" => Vector{Float64}[],
                "lam_kircchoff_active" => Vector{Float64}[],
                "lam_kircchoff_reactive" => Vector{Float64}[],
                "mu_pg_lb" => Vector{Float64}[],
                "mu_pg_ub" => Vector{Float64}[],
                "mu_qg_lb" => Vector{Float64}[],
                "mu_qg_ub" => Vector{Float64}[],
                "mu_sm_fr" => Vector{Float64}[],
                "mu_sm_to" => Vector{Float64}[],
                "lam_ohm_active_fr" => Vector{Float64}[],
                "lam_ohm_active_to" => Vector{Float64}[],
                "lam_ohm_reactive_fr" => Vector{Float64}[],
                "lam_ohm_reactive_to" => Vector{Float64}[],
                "mu_va_diff" => Vector{Float64}[],
            ),
        ),
    )

    return D
end

function add_datapoint!(D, d)
    D_input::Dict{String,Any}    = D["input"]
    D_solution::Dict{String,Any} = D["solution"]

    meta::Dict{String,Any} = d["meta"]
    data::Dict{String,Any} = d["data"]
    res::Dict{String,Any}  = d["res"]
    sol::Dict{String,Any}  = res["solution"]

    N = length(data["bus"])
    E = length(data["branch"])
    L = length(data["load"])
    G = length(data["gen"])

    # TODO: sanity checks

    # Meta info
    seed::Int = meta["seed"]
    push!(D["meta"]["seed"], seed)

    # Input data
    pd::Vector{Float64} = [data["load"]["$l"]["pd"] for l in 1:L]
    qd::Vector{Float64} = [data["load"]["$l"]["qd"] for l in 1:L]
    br::Vector{Bool}    = [Bool(data["branch"]["$e"]["br_status"]) for e in 1:E]
    push!(D_input["pd"], pd)
    push!(D_input["qd"], qd)
    push!(D_input["br_status"], br)

    # Solution meta data
    push!(D["solution"]["meta"]["termination_status"], res["termination_status"])
    push!(D["solution"]["meta"]["primal_status"], res["primal_status"])
    push!(D["solution"]["meta"]["dual_status"], res["dual_status"])
    push!(D["solution"]["meta"]["solve_time"], res["solve_time"])

    # Primal
    vm                     = zeros(Float64, N)
    va                     = zeros(Float64, N)
    pg                     = zeros(Float64, G)
    qg                     = zeros(Float64, G)
    pf                     = zeros(Float64, E)
    qf                     = zeros(Float64, E)
    pt                     = zeros(Float64, E)
    qt                     = zeros(Float64, E)
    # Dual
    mu_vm_lb               = zeros(Float64, N)
    mu_vm_ub               = zeros(Float64, N)
    lam_kircchoff_active   = zeros(Float64, N)
    lam_kircchoff_reactive = zeros(Float64, N)
    mu_pg_lb               = zeros(Float64, G)
    mu_pg_ub               = zeros(Float64, G)
    mu_qg_lb               = zeros(Float64, G)
    mu_qg_ub               = zeros(Float64, G)
    mu_sm_fr               = zeros(Float64, E)
    mu_sm_to               = zeros(Float64, E)
    lam_ohm_active_fr      = zeros(Float64, E)
    lam_ohm_active_to      = zeros(Float64, E)
    lam_ohm_reactive_fr    = zeros(Float64, E)
    lam_ohm_reactive_to    = zeros(Float64, E)
    mu_va_diff             = zeros(Float64, E)

    # extract from solution
    for i in 1:N
        bsol = sol["bus"]["$i"]

        vm[i] = bsol["vm"]
        va[i] = bsol["va"]

        mu_vm_lb[i] = bsol["mu_vm_lb"]
        mu_vm_ub[i] = bsol["mu_vm_ub"]
        lam_kircchoff_active[i] = bsol["lam_pb_active"]
        lam_kircchoff_reactive[i] = bsol["lam_pb_reactive"]
    end
    for g in 1:G
        gsol = sol["gen"]["$g"]

        pg[g] = gsol["pg"]
        qg[g] = gsol["qg"]

        mu_pg_lb[g] = gsol["mu_pg_lb"]
        mu_pg_ub[g] = gsol["mu_pg_ub"]
        mu_qg_lb[g] = gsol["mu_qg_lb"]
        mu_qg_ub[g] = gsol["mu_qg_ub"]
    end
    for e in 1:E
        bsol = sol["branch"]["$e"]

        pf[e] = bsol["pf"]
        qf[e] = bsol["qf"]
        pt[e] = bsol["pt"]
        qt[e] = bsol["qt"]
        
        mu_sm_fr[e] = bsol["mu_sm_fr"]
        mu_sm_to[e] = bsol["mu_sm_to"]
        lam_ohm_active_fr[e] = bsol["lam_ohm_active_fr"]
        lam_ohm_active_to[e] = bsol["lam_ohm_active_to"]
        lam_ohm_reactive_fr[e] = bsol["lam_ohm_reactive_fr"]
        lam_ohm_reactive_to[e] = bsol["lam_ohm_reactive_to"]
        mu_va_diff[e] = bsol["mu_va_diff"]
    end

    # add to Dataset
    push!(D_solution["primal"]["vm"], vm)
    push!(D_solution["primal"]["va"], va)
    push!(D_solution["primal"]["pg"], pg)
    push!(D_solution["primal"]["qg"], qg)
    push!(D_solution["primal"]["pf"], pf)
    push!(D_solution["primal"]["qf"], qf)
    push!(D_solution["primal"]["pt"], pt)
    push!(D_solution["primal"]["qt"], qt)
    push!(D_solution["dual"]["mu_vm_lb"], mu_vm_lb)
    push!(D_solution["dual"]["mu_vm_ub"], mu_vm_ub)
    push!(D_solution["dual"]["lam_kircchoff_active"], lam_kircchoff_active)
    push!(D_solution["dual"]["lam_kircchoff_reactive"], lam_kircchoff_reactive)
    push!(D_solution["dual"]["mu_pg_lb"], mu_pg_lb)
    push!(D_solution["dual"]["mu_pg_ub"], mu_pg_ub)
    push!(D_solution["dual"]["mu_qg_lb"], mu_qg_lb)
    push!(D_solution["dual"]["mu_qg_ub"], mu_qg_ub)
    push!(D_solution["dual"]["mu_sm_fr"], mu_sm_fr)
    push!(D_solution["dual"]["mu_sm_to"], mu_sm_to)
    push!(D_solution["dual"]["lam_ohm_active_fr"], lam_ohm_active_fr)
    push!(D_solution["dual"]["lam_ohm_active_to"], lam_ohm_active_to)
    push!(D_solution["dual"]["lam_ohm_reactive_fr"], lam_ohm_reactive_fr)
    push!(D_solution["dual"]["lam_ohm_reactive_to"], lam_ohm_reactive_to)
    push!(D_solution["dual"]["mu_va_diff"], mu_va_diff)

    return D
end

function save_h5(filename, D)

    h5open(filename, "w") do file
        meta = create_group(file, "meta")
        meta["ref"] = D["meta"]["ref"]
        meta["seed"] = D["meta"]["seed"]

        # Input data
        dat  = create_group(file, "input")
        dat["pd"] = D["input"]["pd"]
        dat["qd"] = D["input"]["qd"]
        dat["br_status"] = D["input"]["br_status"]
        
        # Solution data
        sol  = create_group(file, "solution")
        sol_meta = create_group(sol, "meta")
        sol_meta["termination_status"] = D["solution"]["meta"]["termination_status"]
        sol_meta["primal_status"] = D["solution"]["meta"]["primal_status"]
        sol_meta["dual_status"] = D["solution"]["meta"]["dual_status"]
        sol_meta["solve_time"] = D["solution"]["meta"]["solve_time"]
        
        # The primal/dual solutions are lists of lists, which we convert to matrices for h5
        primal = create_group(sol, "primal")
        primal["vm"] = D["solution"]["primal"]["vm"]
        primal["va"] = D["solution"]["primal"]["va"]
        primal["pg"] = D["solution"]["primal"]["pg"]
        primal["qg"] = D["solution"]["primal"]["qg"]
        primal["pf"] = D["solution"]["primal"]["pf"]
        primal["qf"] = D["solution"]["primal"]["qf"]
        primal["pt"] = D["solution"]["primal"]["pt"]
        primal["qt"] = D["solution"]["primal"]["qt"]

        dual   = create_group(sol, "dual")
        dual["mu_sm_to"] = D["solution"]["dual"]["mu_sm_to"]
        dual["mu_vm_ub"] = D["solution"]["dual"]["mu_vm_ub"]
        dual["mu_pg_ub"] = D["solution"]["dual"]["mu_pg_ub"]
        dual["lam_kircchoff_active"] = D["solution"]["dual"]["lam_kircchoff_active"]
        dual["lam_ohm_reactive_fr"] = D["solution"]["dual"]["lam_ohm_reactive_fr"]
        dual["mu_qg_ub"] = D["solution"]["dual"]["mu_qg_ub"]
        dual["mu_vm_lb"] = D["solution"]["dual"]["mu_vm_lb"]
        dual["mu_sm_fr"] = D["solution"]["dual"]["mu_sm_fr"]
        dual["lam_ohm_active_to"] = D["solution"]["dual"]["lam_ohm_active_to"]
        dual["mu_pg_lb"] = D["solution"]["dual"]["mu_pg_lb"]
        dual["lam_ohm_active_fr"] = D["solution"]["dual"]["lam_ohm_active_fr"]
        dual["mu_qg_lb"] = D["solution"]["dual"]["mu_qg_lb"]
        dual["mu_va_diff"] = D["solution"]["dual"]["mu_va_diff"]
        dual["lam_kircchoff_reactive"] = D["solution"]["dual"]["lam_kircchoff_reactive"]
        dual["lam_ohm_reactive_to"] = D["solution"]["dual"]["lam_ohm_reactive_to"]
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
    # Input data
    dat  = D["input"]
    dat["pd"] = _vecvec2mat(D["input"]["pd"])
    dat["qd"] = _vecvec2mat(D["input"]["qd"])
    dat["br_status"] = _vecvec2mat(D["input"]["br_status"])
        
    # Solution data
    sol  = D["solution"]
        
    # The primal/dual solutions are lists of lists, which we convert to matrices for h5
    primal = sol["primal"]
    primal["vm"] = _vecvec2mat(D["solution"]["primal"]["vm"])
    primal["va"] = _vecvec2mat(D["solution"]["primal"]["va"])
    primal["pg"] = _vecvec2mat(D["solution"]["primal"]["pg"])
    primal["qg"] = _vecvec2mat(D["solution"]["primal"]["qg"])
    primal["pf"] = _vecvec2mat(D["solution"]["primal"]["pf"])
    primal["qf"] = _vecvec2mat(D["solution"]["primal"]["qf"])
    primal["pt"] = _vecvec2mat(D["solution"]["primal"]["pt"])
    primal["qt"] = _vecvec2mat(D["solution"]["primal"]["qt"])

    dual   = sol["dual"]
    dual["mu_sm_to"] = _vecvec2mat(D["solution"]["dual"]["mu_sm_to"])
    dual["mu_vm_ub"] = _vecvec2mat(D["solution"]["dual"]["mu_vm_ub"])
    dual["mu_pg_ub"] = _vecvec2mat(D["solution"]["dual"]["mu_pg_ub"])
    dual["lam_kircchoff_active"] = _vecvec2mat(D["solution"]["dual"]["lam_kircchoff_active"])
    dual["lam_ohm_reactive_fr"] = _vecvec2mat(D["solution"]["dual"]["lam_ohm_reactive_fr"])
    dual["mu_qg_ub"] = _vecvec2mat(D["solution"]["dual"]["mu_qg_ub"])
    dual["mu_vm_lb"] = _vecvec2mat(D["solution"]["dual"]["mu_vm_lb"])
    dual["mu_sm_fr"] = _vecvec2mat(D["solution"]["dual"]["mu_sm_fr"])
    dual["lam_ohm_active_to"] = _vecvec2mat(D["solution"]["dual"]["lam_ohm_active_to"])
    dual["mu_pg_lb"] = _vecvec2mat(D["solution"]["dual"]["mu_pg_lb"])
    dual["lam_ohm_active_fr"] = _vecvec2mat(D["solution"]["dual"]["lam_ohm_active_fr"])
    dual["mu_qg_lb"] = _vecvec2mat(D["solution"]["dual"]["mu_qg_lb"])
    dual["mu_va_diff"] = _vecvec2mat(D["solution"]["dual"]["mu_va_diff"])
    dual["lam_kircchoff_reactive"] = _vecvec2mat(D["solution"]["dual"]["lam_kircchoff_reactive"])
    dual["lam_ohm_reactive_to"] = _vecvec2mat(D["solution"]["dual"]["lam_ohm_reactive_to"])

    return nothing
end

_vecvec2mat(V) = reduce(hcat, V)

function parse_result_folder(result_folder::String, export_folder::String;
    show_progress::Bool=true,
    batch_size::Int=0,
    force_process::Bool=true,
)

    all_files = filter(s -> endswith(s, r".json|.json.gz"), readdir(result_folder))
    sort!(all_files)
    N = length(all_files)

    # Read first file
    fname = all_files[1]
    ref = load_json(joinpath(result_folder, fname))["meta"]["ref"]
    data = make_basic_network(pglib(ref))

    # If no batch size was provided, process everything in one go
    batch_size = (batch_size <= 0) ? N : batch_size
    F = partition(all_files, batch_size)
    println("Processing $N results files from $(result_folder)")
    println("Chunk size: $(batch_size) ($(length(F)) chunks)")
    println("Results will be exported to $(joinpath(export_folder, ref * "<chunk_index>.h5."))")

    p = Progress(N; enabled=show_progress)

    for (b, files) in enumerate(F)
        f5name = joinpath(export_folder, ref * "_$(b).h5")
        !force_process && isfile(f5name) && continue

        D = initialize_res(data)

        # Process all files by chunks
        L = ReentrantLock()
        num_read_error = Threads.Atomic{Int}(0)
        @threads for fname in files
            next!(p)
            local d = try
                load_json(joinpath(result_folder, fname))
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
    p = sortperm(D["meta"]["seed"])

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

if abspath(PROGRAM_FILE) == @__FILE__
    result_folder = ARGS[1]
    export_folder = ARGS[2]
    batch_size    = parse(Int, ARGS[3])
    parse_result_folder(result_folder, export_folder;
        show_progress=false,
        batch_size,
        force_process=false,
    )
end
