using Random

function check_dataset(D, result_folder; S::Int=0)
    ref = D["meta"]["ref"]
    data0 = make_basic_network(pglib(ref))

    N = length(data0["bus"])
    E = length(data0["branch"])
    G = length(data0["gen"])
    L = length(data0["load"])

    rng = MersenneTwister(42)
    nseeds = length(D["meta"]["seed"])
    S = (S <= 0) ? nseeds : min(S, nseeds)
    J = sort(D["meta"]["seed"][randperm(rng, nseeds)[1:S]])

    for j in J
        s = D["meta"]["seed"][j]
        local d = load_json(joinpath(result_folder, ref * "_s$(s).json.gz"))
        data = d["data"]
        sol  = d["res"]["solution"]

        @assert d["meta"]["seed"] == s
        
        # check inputs
        @assert D["input"]["pd"][:, j] == [d["data"]["load"]["$l"]["pd"] for l in 1:L] s
        @assert D["input"]["qd"][:, j] == [d["data"]["load"]["$l"]["qd"] for l in 1:L] s
        
        # check primal solution
        @assert D["solution"]["primal"]["pg"][:, j] == [sol["gen"]["$g"]["pg"] for g in 1:G] s
        @assert D["solution"]["primal"]["qg"][:, j] == [sol["gen"]["$g"]["qg"] for g in 1:G] s

        # check dual solution
        @assert D["solution"]["dual"]["mu_vm_lb"][:, j] == [sol["bus"]["$i"]["mu_vm_lb"] for i in 1:N] s
        @assert D["solution"]["dual"]["mu_vm_ub"][:, j] == [sol["bus"]["$i"]["mu_vm_ub"] for i in 1:N] s
    end

    @info "Checked $(S) instances and found no issue."

    return nothing
end
