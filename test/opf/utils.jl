function test_arcs_from(data, ref)
    arcs_from = OPFGenerator.make_arcs_from(data)
    @test arcs_from == sort(ref[:arcs_from], by=x->x[1])
end

function test_arcs_to(data, ref)
    arcs_to = OPFGenerator.make_arcs_to(data)
    @test arcs_to == sort(ref[:arcs_to], by=x->x[1])
end

function test_arcs(data, ref)
    arcs = OPFGenerator.make_arcs(data)
    @test length(arcs) == length(ref[:arcs])
    # `arcs` is ordered like [arcs_from; arcs_to]
    # ideally sort ref[:arcs] to match but sorting both suffices
    @test sort(arcs, by=x->x[1]) == sort(ref[:arcs], by=x->(x[1]))
end

function test_bus_loads(data, ref)
    bus_loads = OPFGenerator.make_bus_loads(data)
    rbus_loads = ref[:bus_loads]
    for (k, v) in rbus_loads
        rbus_loads[k] = sort(v)
    end
    @test bus_loads == [rbus_loads[b] for b in 1:length(data["bus"])]
end

function test_bus_shunts(data, ref)
    bus_shunts = OPFGenerator.make_bus_shunts(data)
    rbus_shunts = ref[:bus_shunts]
    for (k, v) in rbus_shunts
        rbus_shunts[k] = sort(v)
    end
    @test bus_shunts == [rbus_shunts[b] for b in 1:length(data["bus"])]
end

function test_bus_gens(data, ref)
    bus_gens = OPFGenerator.make_bus_gens(data)
    rbus_gens = ref[:bus_gens]
    for (k, v) in rbus_gens
        rbus_gens[k] = sort(v)
    end
    @test bus_gens == [rbus_gens[b] for b in 1:length(data["bus"])]
end

function test_bus_arcs(data, ref)
    bus_arcs = OPFGenerator.make_bus_arcs(data)
    rbus_arcs = ref[:bus_arcs]
    for (k, v) in rbus_arcs
        rbus_arcs[k] = sort(v, by=x->x[1])
    end
    @test bus_arcs == [rbus_arcs[b] for b in 1:length(data["bus"])]
end

function test_ref_buses(data, ref)
    ref_buses = OPFGenerator.make_ref_buses(data)
    @test ref_buses == sort(collect(keys(ref[:ref_buses])))
end

function test_active_branches(data, ref)
    active_branches = OPFGenerator.make_active_branches(data)
    @test active_branches == sort(collect(keys(ref[:branch])))

    data_ = deepcopy(data)
    data_["branch"]["1"]["br_status"] = 0
    n1_active_branches = OPFGenerator.make_active_branches(data)
    n1_ref = PM.build_ref(data)[:it][:pm][:nw][0]
    @test n1_active_branches == sort(collect(keys(n1_ref[:branch])))
end

function test_active_generators(data, ref)
    active_generators = OPFGenerator.make_active_generators(data)
    @test active_generators == sort(collect(keys(ref[:gen])))

    data_ = deepcopy(data)
    data_["gen"]["1"]["gen_status"] = 0
    n1_active_generators = OPFGenerator.make_active_generators(data)
    n1_ref = PM.build_ref(data)[:it][:pm][:nw][0]
    @test n1_active_generators == sort(collect(keys(n1_ref[:gen])))
end

# function test_fix_variables(data, ref)
#     # TODO
# end

function test_calc_buspair_voltage_bounds(data, ref)
    data_ = deepcopy(data)
    calc_buspair_voltage_bounds!(data_)

    wr_min = Dict((i,j) => data_["buspair"][(i,j)]["wr_min"] for (i,j) in keys(data_["buspair"]))
    wr_max = Dict((i,j) => data_["buspair"][(i,j)]["wr_max"] for (i,j) in keys(data_["buspair"]))
    wi_min = Dict((i,j) => data_["buspair"][(i,j)]["wi_min"] for (i,j) in keys(data_["buspair"]))
    wi_max = Dict((i,j) => data_["buspair"][(i,j)]["wi_max"] for (i,j) in keys(data_["buspair"]))

    buspairs = PowerModels.calc_buspair_parameters(data["bus"], data["branch"])
    wr_min_pm, wr_max_pm, wi_min_pm, wi_max_pm = PowerModels.calc_buspair_voltage_bounds(buspairs)

    @test wr_min == wr_min_pm
    @test wr_max == wr_max_pm
    @test wi_min == wi_min_pm
    @test wi_max == wi_max_pm
end

function test_utils(casename)
    data = PM.make_basic_network(pglib(casename))
    rref = PM.build_ref(data)[:it][:pm][:nw][0]
    @testset test_arcs_from(data, rref)
    @testset test_arcs_to(data, rref)
    @testset test_arcs(data, rref)
    @testset test_bus_loads(data, rref)
    @testset test_bus_shunts(data, rref)
    @testset test_bus_gens(data, rref)
    @testset test_bus_arcs(data, rref)
    @testset test_ref_buses(data, rref)
    @testset test_active_branches(data, rref)
    @testset test_active_generators(data, rref)
end

@testset "Utils" begin
    @testset "$(casename)" for casename in PGLIB_CASES
        test_utils(casename)
    end
end