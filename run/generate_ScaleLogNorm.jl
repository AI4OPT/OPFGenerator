using ACOPFGenerator, 
        PowerModels, 
        StableRNGs, 
        PGLib, 
        JuMP, 
        Ipopt, 
        Distributions

"""
    create_seed_range(seed, total[, max])

Divides a large range over many computing nodes 
by generating a range of int seeds in [1, `total`],
with size proportional to `total`/`max`, where 
`max` is the maximum nodes allowed per job on a 
computing cluster, and `total` is the amount of 
separate jobs desired.

Example:
create_seed_range(1, 1000, 500) 
>  1:2

create_seed_range(403, 50000, 500)
>  40201:40300

create_seed_range(210, 400, 500)
>  210:210
"""
function create_seed_range(seed, total, max = 500)
    if total > max
        remainder = total % max
        cap = max - remainder
        scale = total ÷ max
        base = cap*scale

        if seed > cap
            scale = scale+1
            upper = (x-cap)*(scale) + base
        else
            upper = (scale*seed)
        end 

        lower = upper - scale + 1
        return range(lower, upper)
    else 
        return range(seed, seed)
    end
end

function main(args)
    # name of PGLib instance
    pglib_file = args[1]

    # lower bound of base load factor range
    l = parse(Float64, args[2])

    # upper bound of base load factor range
    u = parse(Float64, args[3])

    # noise level
    v = parse(Float64, args[4])

    # random seed passed in as an integer
    seed = parse(Int64, args[5])

    # where generated .json will be saved
    save_loc = args[6]

    # number of augmented instances
    total = parse(Int64, args[7])

    # cluster limit on parallel jobs  
    max = parse(Int64, args[8])

    seed_range = create_seed_range(seed, total, max)
    
    data = make_basic_network(pglib(pglib_file))
    sampler = ScaleLogNorm(data, l, u, v);
    opf = SimpleOPFSampler(data, sampler);
        
    range_save(seed_range, opf, save_loc); 
end 

main(ARGS)