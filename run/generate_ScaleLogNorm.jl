using Pkg;
Pkg.activate(".")
using ACOPFGenerator, PowerModels, StableRNGs, PGLib, JuMP, Ipopt, Distributions

function main(args)
    pglib_file = args[1]
    l = parse(Float64, args[2])
    u = parse(Float64, args[3])
    v = parse(Float64, args[4])
    num_loads = parse(Int64, args[5])
    iter = parse(Int64, args[6])
    save_file = args[7]
    
    data = make_basic_network(pglib(pglib_file))
    sampler = ScaleLogNorm(data, l, u, v);
    opf = SimpleOPFSampler(data, sampler);
    
    range_save(range(iter, iter), opf, save_file);
    
end 

main(ARGS)