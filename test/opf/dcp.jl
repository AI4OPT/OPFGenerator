_test_dcopf(casename::String) = _test_dcopf(PM.make_basic_network(pglib(casename)))

"""
    _test_dcopf(data)

Solve DCOPF problem and compare to PowerModels' solution.
"""
function _test_dcopf(data::Dict)
    data["basic_network"] || error("Input data must be in basic format to test")
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])

    solver = optimizer_with_attributes(Ipopt.Optimizer,
        "print_level" => 0,
    )
    res_pm = solve_dc_opf(data, solver)
    sol_pm = res_pm["solution"]

    # build and solve DCOPF 
    dcopf = OPFGenerator.build_opf(PM.DCPPowerModel, data, solver)
    model = dcopf.model
    optimize!(model)
    res = OPFGenerator.extract_result(dcopf)
    @test res["opf_model"] == "DCPPowerModel"
    sol = res["solution"]

    # Check that we get consistent results with PowerModels
    @test res["termination_status"] == LOCALLY_SOLVED
    @test isapprox(res["objective"], res_pm["objective"], atol=1e-6, rtol=1e-6)

    # Force PM solution into our model, and check that the solution is feasible
    # TODO: use JuMP.primal_feasibility_report instead
    #    (would require extracting a variable => value Dict)
    var2val_pm = Dict(
        :pg => Float64[sol_pm["gen"]["$g"]["pg"] for g in 1:G],
        :va => Float64[sol_pm["bus"]["$i"]["va"] for i in 1:N],
        :pf => Float64[sol_pm["branch"]["$e"]["pf"] for e in 1:E],
    )
    # reorder pf according to branch order in our model
    dcopf_branch_order = [key[1][1] for key in keys(model[:pf])][1:E]
    var2val_pm[:pf] = var2val_pm[:pf][dcopf_branch_order]
    var2val_pm[:pf] = vcat(var2val_pm[:pf], -var2val_pm[:pf])

    @constraint(model, var2val_pm[:pg] .<= model[:pg] .<= var2val_pm[:pg])
    @constraint(model, var2val_pm[:va] .<= model[:va] .<= var2val_pm[:va])
    @constraint(model, var2val_pm[:pf] .<= model[:pf] .<= var2val_pm[:pf])

    optimize!(model)
    @test termination_status(model) ∈ [LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    @test primal_status(model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]

    return nothing
end

@testset "DCOPF" begin
    @testset "$casename" for casename in ["14_ieee", "30_ieee", "57_ieee", "89_pegase", "118_ieee"]
        _test_dcopf("pglib_opf_case$(casename)")
    end
end

