"""
    bundle_solve(i::Int64, sampler::SimpleOPFSampler)

Augments data stored in `sampler` using an RNG with seed `i`. Returns a dictionary of the 
augmented data, the instance solved with PowerModels.ac_opf_solution(), and metadata for the solution. 
"""
function bundle_solve(i::Int64, sampler::SimpleOPFSampler)
    rng = StableRNG(i)
    agmtd_data = sample(rng, sampler)
    sol = solve_ac_opf(agmtd_data, Ipopt.Optimizer)
    meta = Dict("network" => sampler.data["name"],
                "seed" => Dict(i => rng),
                "augment_method" => Dict(typeof(sampler.load_sampler) => sampler.load_sampler))
    return Dict("agmtd_data" => agmtd_data, "sol" => sol, "meta" => meta)
end 


"""
    range_save(iters::AbstractRange, sampler::SimpleOPFSampler, loc::AbstractString)

Each value in range `iters` saves one json.gz file to location `loc` containing the output of
`simple_solve`'s dictionary with inputted data contained in `sampler`. 
"""
function range_save(iters::AbstractRange, sampler::SimpleOPFSampler, loc::AbstractString)
    isdir(loc) || error("folder location does not exist. '/' may be needed at beginning of save location")

    for i in iters
        bundle = bundle_solve(i, sampler)
        name_long = bundle["meta"]["network"]
        name_short = net_name[15:length(net_name)]
        fn = name_short * "_" * string(i) * ".json.gz"
        save_path = loc * "/" * fn
        save_json(save_path, bundle)
    end 
end