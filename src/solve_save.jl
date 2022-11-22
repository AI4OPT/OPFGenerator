"""
    augment_solve(seed::AbstractRNG, sampler::SimpleOPFSampler)

Augments data stored in `sampler` using an RNG `rng`. Returns a dictionary of the 
augmented data, the instance solved extracted with acopf._extract_solution(), and 
metadata for the solution. 
"""
function augment_solve(rng::AbstractRNG, sampler::SimpleOPFSampler)
    data = sample(rng, sampler)
    model = build_acopf(data, Ipopt.Optimizer)
    optimize!(model)
    sol = _extract_solution(model, data)
    meta = Dict(
        "network" => sampler.data["name"],
        "seed" => rng,
        "augment_method" => Dict(
            "name" => "$(typeof(sampler.load_sampler))",
            "base_load_factors" => (sampler.load_sampler.d0.a, sampler.load_sampler.d0.b),
            "noise_level" => sampler.load_sampler.ds[1].σ
        )
    )
    return Dict("data" => data, "sol" => sol, "meta" => meta)
end 


"""
    range_save(iters::AbstractRange, sampler::SimpleOPFSampler, loc::AbstractString)

Each value in range `iters` saves one json.gz file to location `loc` containing the output of
`simple_solve`'s dictionary with inputted data contained in `sampler`. 
"""
function range_save(iters::AbstractRange, sampler::SimpleOPFSampler, loc::AbstractString)
    isdir(loc) || error("folder location does not exist. '/' may be needed at beginning of save location")

    for seed in iters
        rng = StableRNG(seed)
        bundle = augment_solve(rng, sampler)
        net_name = bundle["meta"]["network"]
        fn = net_name * "_" * string(seed) * ".json.gz"
        save_path = joinpath(loc, fn)
        save_json(save_path, bundle)
    end 
end
