using Pkg;
Pkg.activate(".")
using ACOPFGenerator, PowerModels, StableRNGs, PGLib, JuMP, Ipopt, Distributions

function create_seed_range(x, total, max = 500)
    if total > max
        remainder = total % max
        cap = max - remainder
        scale = total ÷ max
        base = cap*scale

        if x > cap
            scale = scale+1
            upper = (x-cap)*(scale) + base
        else
            upper = (scale*x)
        end 

        lower = upper - scale + 1
        return range(lower, upper)
    else 
        return range(x, x)
    end
end

function main(args)

    pglib_file = args[1]
    l = parse(Float64, args[2])
    u = parse(Float64, args[3])
    v = parse(Float64, args[4])
    seed = parse(Int64, args[5])
    save_file = args[6]
    total = parse(Int64, args[7])
    max = parse(Int64, args[8])

    seed_range = create_seed_range(seed, total, max)
    
    data = make_basic_network(pglib(pglib_file))
    sampler = ScaleLogNorm(data, l, u, v);
    opf = SimpleOPFSampler(data, sampler);
        
    range_save(seed_range, opf, save_file); 
end 

main(ARGS)