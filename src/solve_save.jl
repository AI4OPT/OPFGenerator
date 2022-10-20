"""
    bundle_solve(i::Int64, sampler::SimpleOPFSampler)

Augments data stored in `sampler` using an RNG with seed `i`. Returns a dictionary of the 
augmented data, the instance solved with PowerModels.ac_opf_solution(), and metadata for the solution. 
"""
function bundle_solve(seed::Int64, sampler::SimpleOPFSampler)
    rng = StableRNG(seed)
    agmtd_data = sample(rng, sampler)
    sol = solve_ac_opf(agmtd_data, Ipopt.Optimizer)
    meta = Dict(
        "network" => sampler.data["name"],
        "seed" => Dict(seed => rng),
        "augment_method" => Dict(
            "name" => "$(typeof(sampler.load_sampler))",
            "base_load_factors" => (sampler.load_sampler.d0.a, sampler.load_sampler.d0.b),
            "noise_level" => sampler.load_sampler.ds[1].σ
        )
    )
    return Dict("agmtd_data" => agmtd_data, "sol" => sol, "meta" => meta)
end 


"""
    range_save(iters::AbstractRange, sampler::SimpleOPFSampler, loc::AbstractString)

Each value in range `iters` saves one json.gz file to location `loc` containing the output of
`simple_solve`'s dictionary with inputted data contained in `sampler`. 
"""
function range_save(iters::AbstractRange, sampler::SimpleOPFSampler, loc::AbstractString)
    isdir(loc) || error("folder location does not exist. '/' may be needed at beginning of save location")

    for seed in iters
        bundle = bundle_solve(seed, sampler)
        net_name = bundle["meta"]["network"]
        name_short = net_name[15:length(net_name)]
        fn = name_short * "_" * string(seed) * ".json.gz"
        save_path = loc * "/" * fn
        save_json(save_path, bundle)
    end 
end