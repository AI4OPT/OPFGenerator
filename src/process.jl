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
        "acopf_solution" => Dict{String,Any}(
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
        "dcopf_solution" => Dict{String,Any}(
            "meta" => Dict{String,Any}(
                "termination_status" => String[],
                "primal_status" => String[],
                "dual_status" => String[],
                "solve_time" => Float64[],
            ),
            "primal" => Dict{String,Any}(
                "va" => Vector{Float64}[],
                "pg" => Vector{Float64}[],
                "pf" => Vector{Float64}[],
            ),
            "dual" => Dict{String,Any}(
                "lam_kirchhoff" => Vector{Float64}[],
                "mu_pg_lb" => Vector{Float64}[],
                "mu_pg_ub" => Vector{Float64}[],
                "lam_ohm" => Vector{Float64}[],
                "mu_va_diff" => Vector{Float64}[],
            ),
        ),
    )

    return D
end

function add_datapoint!(D, d)
    D_input::Dict{String,Any}    = D["input"]
    meta::Dict{String,Any} = d["meta"]
    data::Dict{String,Any} = d["data"]

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

    # Add ACOPF Solution
    D_acopf_solution::Dict{String,Any} = D["acopf_solution"]
    acopf_res::Dict{String,Any}  = d["acopf_res"]
    acopf_sol::Dict{String,Any}  = acopf_res["solution"]

    # ACOPF Solution meta data
    push!(D_acopf_solution["meta"]["termination_status"], acopf_res["termination_status"])
    push!(D_acopf_solution["meta"]["primal_status"], acopf_res["primal_status"])
    push!(D_acopf_solution["meta"]["dual_status"], acopf_res["dual_status"])
    push!(D_acopf_solution["meta"]["solve_time"], acopf_res["solve_time"])

    # ACOPF Primal
    vm                     = zeros(Float64, N)
    va                     = zeros(Float64, N)
    pg                     = zeros(Float64, G)
    qg                     = zeros(Float64, G)
    pf                     = zeros(Float64, E)
    qf                     = zeros(Float64, E)
    pt                     = zeros(Float64, E)
    qt                     = zeros(Float64, E)
    # ACOPF Dual
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

    # extract from ACOPF solution
    for i in 1:N
        bsol = acopf_sol["bus"]["$i"]

        vm[i] = bsol["vm"]
        va[i] = bsol["va"]

        mu_vm_lb[i] = bsol["mu_vm_lb"]
        mu_vm_ub[i] = bsol["mu_vm_ub"]
        lam_kircchoff_active[i] = bsol["lam_pb_active"]
        lam_kircchoff_reactive[i] = bsol["lam_pb_reactive"]
    end
    for g in 1:G
        gsol = acopf_sol["gen"]["$g"]

        pg[g] = gsol["pg"]
        qg[g] = gsol["qg"]

        mu_pg_lb[g] = gsol["mu_pg_lb"]
        mu_pg_ub[g] = gsol["mu_pg_ub"]
        mu_qg_lb[g] = gsol["mu_qg_lb"]
        mu_qg_ub[g] = gsol["mu_qg_ub"]
    end
    for e in 1:E
        bsol = acopf_sol["branch"]["$e"]

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

    # add ACOPF to Dataset
    push!(D_acopf_solution["primal"]["vm"], vm)
    push!(D_acopf_solution["primal"]["va"], va)
    push!(D_acopf_solution["primal"]["pg"], pg)
    push!(D_acopf_solution["primal"]["qg"], qg)
    push!(D_acopf_solution["primal"]["pf"], pf)
    push!(D_acopf_solution["primal"]["qf"], qf)
    push!(D_acopf_solution["primal"]["pt"], pt)
    push!(D_acopf_solution["primal"]["qt"], qt)
    push!(D_acopf_solution["dual"]["mu_vm_lb"], mu_vm_lb)
    push!(D_acopf_solution["dual"]["mu_vm_ub"], mu_vm_ub)
    push!(D_acopf_solution["dual"]["lam_kircchoff_active"], lam_kircchoff_active)
    push!(D_acopf_solution["dual"]["lam_kircchoff_reactive"], lam_kircchoff_reactive)
    push!(D_acopf_solution["dual"]["mu_pg_lb"], mu_pg_lb)
    push!(D_acopf_solution["dual"]["mu_pg_ub"], mu_pg_ub)
    push!(D_acopf_solution["dual"]["mu_qg_lb"], mu_qg_lb)
    push!(D_acopf_solution["dual"]["mu_qg_ub"], mu_qg_ub)
    push!(D_acopf_solution["dual"]["mu_sm_fr"], mu_sm_fr)
    push!(D_acopf_solution["dual"]["mu_sm_to"], mu_sm_to)
    push!(D_acopf_solution["dual"]["lam_ohm_active_fr"], lam_ohm_active_fr)
    push!(D_acopf_solution["dual"]["lam_ohm_active_to"], lam_ohm_active_to)
    push!(D_acopf_solution["dual"]["lam_ohm_reactive_fr"], lam_ohm_reactive_fr)
    push!(D_acopf_solution["dual"]["lam_ohm_reactive_to"], lam_ohm_reactive_to)
    push!(D_acopf_solution["dual"]["mu_va_diff"], mu_va_diff)

    # Add DCOPF solution
    D_dcopf_solution::Dict{String,Any} = D["dcopf_solution"]
    dcopf_res::Dict{String,Any}  = d["dcopf_res"]
    dcopf_sol::Dict{String,Any}  = dcopf_res["solution"]

    # DCOPF solution meta data
    push!(D_dcopf_solution["meta"]["termination_status"], dcopf_res["termination_status"])
    push!(D_dcopf_solution["meta"]["primal_status"], dcopf_res["primal_status"])
    push!(D_dcopf_solution["meta"]["dual_status"], dcopf_res["dual_status"])
    push!(D_dcopf_solution["meta"]["solve_time"], dcopf_res["solve_time"])

    # DCOPF Primal
    va            = zeros(Float64, N)
    pg            = zeros(Float64, G)
    pf            = zeros(Float64, E)
    # DCOPF Dual
    lam_kirchhoff = zeros(Float64, N)
    mu_pg_lb      = zeros(Float64, G)
    mu_pg_ub      = zeros(Float64, G)
    lam_ohm       = zeros(Float64, E)
    mu_va_diff    = zeros(Float64, E)

    # extract from DCOPF solution
    for i in 1:N
        bsol = dcopf_sol["bus"]["$i"]

        va[i] = bsol["va"]

        lam_kirchhoff[i] = bsol["lam_pb"]
    end
    for g in 1:G
        gsol = dcopf_sol["gen"]["$g"]

        pg[g] = gsol["pg"]

        mu_pg_lb[g] = gsol["mu_pg_lb"]
        mu_pg_ub[g] = gsol["mu_pg_ub"]
    end
    for e in 1:E
        bsol = dcopf_sol["branch"]["$e"]

        pf[e] = bsol["pf"]

        lam_ohm[e] = bsol["lam_ohm"]
        mu_va_diff[e] = bsol["mu_va_diff"]
    end

    # DCOPF add to Dataset
    push!(D_dcopf_solution["primal"]["va"], va)
    push!(D_dcopf_solution["primal"]["pg"], pg)
    push!(D_dcopf_solution["primal"]["pf"], pf)
    push!(D_dcopf_solution["dual"]["lam_kirchhoff"], lam_kirchhoff)
    push!(D_dcopf_solution["dual"]["mu_pg_lb"], mu_pg_lb)
    push!(D_dcopf_solution["dual"]["mu_pg_ub"], mu_pg_ub)
    push!(D_dcopf_solution["dual"]["lam_ohm"], lam_ohm)
    push!(D_dcopf_solution["dual"]["mu_va_diff"], mu_va_diff)

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
        
        # ACOPF Solution data
        acopf_sol  = create_group(file, "acopf_solution")
        acopf_sol_meta = create_group(acopf_sol, "meta")
        acopf_sol_meta["termination_status"] = D["acopf_solution"]["meta"]["termination_status"]
        acopf_sol_meta["primal_status"] = D["acopf_solution"]["meta"]["primal_status"]
        acopf_sol_meta["dual_status"] = D["acopf_solution"]["meta"]["dual_status"]
        acopf_sol_meta["solve_time"] = D["acopf_solution"]["meta"]["solve_time"]
        
        # The primal/dual solutions are lists of lists, which we convert to matrices for h5
        acopf_primal = create_group(acopf_sol, "primal")
        acopf_primal["vm"] = D["acopf_solution"]["primal"]["vm"]
        acopf_primal["va"] = D["acopf_solution"]["primal"]["va"]
        acopf_primal["pg"] = D["acopf_solution"]["primal"]["pg"]
        acopf_primal["qg"] = D["acopf_solution"]["primal"]["qg"]
        acopf_primal["pf"] = D["acopf_solution"]["primal"]["pf"]
        acopf_primal["qf"] = D["acopf_solution"]["primal"]["qf"]
        acopf_primal["pt"] = D["acopf_solution"]["primal"]["pt"]
        acopf_primal["qt"] = D["acopf_solution"]["primal"]["qt"]

        acopf_dual   = create_group(acopf_sol, "dual")
        acopf_dual["mu_sm_to"] = D["acopf_solution"]["dual"]["mu_sm_to"]
        acopf_dual["mu_vm_ub"] = D["acopf_solution"]["dual"]["mu_vm_ub"]
        acopf_dual["mu_pg_ub"] = D["acopf_solution"]["dual"]["mu_pg_ub"]
        acopf_dual["lam_kircchoff_active"] = D["acopf_solution"]["dual"]["lam_kircchoff_active"]
        acopf_dual["lam_ohm_reactive_fr"] = D["acopf_solution"]["dual"]["lam_ohm_reactive_fr"]
        acopf_dual["mu_qg_ub"] = D["acopf_solution"]["dual"]["mu_qg_ub"]
        acopf_dual["mu_vm_lb"] = D["acopf_solution"]["dual"]["mu_vm_lb"]
        acopf_dual["mu_sm_fr"] = D["acopf_solution"]["dual"]["mu_sm_fr"]
        acopf_dual["lam_ohm_active_to"] = D["acopf_solution"]["dual"]["lam_ohm_active_to"]
        acopf_dual["mu_pg_lb"] = D["acopf_solution"]["dual"]["mu_pg_lb"]
        acopf_dual["lam_ohm_active_fr"] = D["acopf_solution"]["dual"]["lam_ohm_active_fr"]
        acopf_dual["mu_qg_lb"] = D["acopf_solution"]["dual"]["mu_qg_lb"]
        acopf_dual["mu_va_diff"] = D["acopf_solution"]["dual"]["mu_va_diff"]
        acopf_dual["lam_kircchoff_reactive"] = D["acopf_solution"]["dual"]["lam_kircchoff_reactive"]
        acopf_dual["lam_ohm_reactive_to"] = D["acopf_solution"]["dual"]["lam_ohm_reactive_to"]

        # DCOPF Solution data
        dcopf_sol  = create_group(file, "dcopf_solution")
        dcopf_sol_meta = create_group(dcopf_sol, "meta")
        dcopf_sol_meta["termination_status"] = D["dcopf_solution"]["meta"]["termination_status"]
        dcopf_sol_meta["primal_status"] = D["dcopf_solution"]["meta"]["primal_status"]
        dcopf_sol_meta["dual_status"] = D["dcopf_solution"]["meta"]["dual_status"]
        dcopf_sol_meta["solve_time"] = D["dcopf_solution"]["meta"]["solve_time"]
        
        # The primal/dual solutions are lists of lists, which we convert to matrices for h5
        dcopf_primal = create_group(dcopf_sol, "primal")
        dcopf_primal["va"] = D["dcopf_solution"]["primal"]["va"]
        dcopf_primal["pf"] = D["dcopf_solution"]["primal"]["pf"]
        dcopf_primal["pg"] = D["dcopf_solution"]["primal"]["pg"]

        dcopf_dual   = create_group(dcopf_sol, "dual")
        dcopf_dual["lam_kirchhoff"] = D["dcopf_solution"]["dual"]["lam_kirchhoff"]
        dcopf_dual["mu_va_diff"] = D["dcopf_solution"]["dual"]["mu_va_diff"]
        dcopf_dual["lam_ohm"] = D["dcopf_solution"]["dual"]["lam_ohm"]
        dcopf_dual["mu_pg_lb"] = D["dcopf_solution"]["dual"]["mu_pg_lb"]
        dcopf_dual["mu_pg_ub"] = D["dcopf_solution"]["dual"]["mu_pg_ub"]
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
        
    # ACOPF Solution data
    acopf_sol  = D["acopf_solution"]
        
    # The primal/dual solutions are lists of lists, which we convert to matrices for h5
    acopf_primal = acopf_sol["primal"]
    acopf_primal["vm"] = _vecvec2mat(D["acopf_solution"]["primal"]["vm"])
    acopf_primal["va"] = _vecvec2mat(D["acopf_solution"]["primal"]["va"])
    acopf_primal["pg"] = _vecvec2mat(D["acopf_solution"]["primal"]["pg"])
    acopf_primal["qg"] = _vecvec2mat(D["acopf_solution"]["primal"]["qg"])
    acopf_primal["pf"] = _vecvec2mat(D["acopf_solution"]["primal"]["pf"])
    acopf_primal["qf"] = _vecvec2mat(D["acopf_solution"]["primal"]["qf"])
    acopf_primal["pt"] = _vecvec2mat(D["acopf_solution"]["primal"]["pt"])
    acopf_primal["qt"] = _vecvec2mat(D["acopf_solution"]["primal"]["qt"])

    acopf_dual   = acopf_sol["dual"]
    acopf_dual["mu_sm_to"] = _vecvec2mat(D["acopf_solution"]["dual"]["mu_sm_to"])
    acopf_dual["mu_vm_ub"] = _vecvec2mat(D["acopf_solution"]["dual"]["mu_vm_ub"])
    acopf_dual["mu_pg_ub"] = _vecvec2mat(D["acopf_solution"]["dual"]["mu_pg_ub"])
    acopf_dual["lam_kircchoff_active"] = _vecvec2mat(D["acopf_solution"]["dual"]["lam_kircchoff_active"])
    acopf_dual["lam_ohm_reactive_fr"] = _vecvec2mat(D["acopf_solution"]["dual"]["lam_ohm_reactive_fr"])
    acopf_dual["mu_qg_ub"] = _vecvec2mat(D["acopf_solution"]["dual"]["mu_qg_ub"])
    acopf_dual["mu_vm_lb"] = _vecvec2mat(D["acopf_solution"]["dual"]["mu_vm_lb"])
    acopf_dual["mu_sm_fr"] = _vecvec2mat(D["acopf_solution"]["dual"]["mu_sm_fr"])
    acopf_dual["lam_ohm_active_to"] = _vecvec2mat(D["acopf_solution"]["dual"]["lam_ohm_active_to"])
    acopf_dual["mu_pg_lb"] = _vecvec2mat(D["acopf_solution"]["dual"]["mu_pg_lb"])
    acopf_dual["lam_ohm_active_fr"] = _vecvec2mat(D["acopf_solution"]["dual"]["lam_ohm_active_fr"])
    acopf_dual["mu_qg_lb"] = _vecvec2mat(D["acopf_solution"]["dual"]["mu_qg_lb"])
    acopf_dual["mu_va_diff"] = _vecvec2mat(D["acopf_solution"]["dual"]["mu_va_diff"])
    acopf_dual["lam_kircchoff_reactive"] = _vecvec2mat(D["acopf_solution"]["dual"]["lam_kircchoff_reactive"])
    acopf_dual["lam_ohm_reactive_to"] = _vecvec2mat(D["acopf_solution"]["dual"]["lam_ohm_reactive_to"])

    # DCOPF Solution data
    dcopf_sol  = D["dcopf_solution"]
    
    # The primal/dual solutions are lists of lists, which we convert to matrices for h5
    dcopf_primal = dcopf_sol["primal"]
    dcopf_primal["va"] = _vecvec2mat(D["dcopf_solution"]["primal"]["va"])
    dcopf_primal["pg"] = _vecvec2mat(D["dcopf_solution"]["primal"]["pg"])
    dcopf_primal["pf"] = _vecvec2mat(D["dcopf_solution"]["primal"]["pf"])

    dcopf_dual   = dcopf_sol["dual"]
    dcopf_dual["lam_kirchhoff"] = _vecvec2mat(D["dcopf_solution"]["dual"]["lam_kirchhoff"])
    dcopf_dual["mu_va_diff"] = _vecvec2mat(D["dcopf_solution"]["dual"]["mu_va_diff"])
    dcopf_dual["lam_ohm"] = _vecvec2mat(D["dcopf_solution"]["dual"]["lam_ohm"])
    dcopf_dual["mu_pg_lb"] = _vecvec2mat(D["dcopf_solution"]["dual"]["mu_pg_lb"])
    dcopf_dual["mu_pg_ub"] = _vecvec2mat(D["dcopf_solution"]["dual"]["mu_pg_ub"])

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
