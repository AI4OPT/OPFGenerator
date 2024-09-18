using OPFGenerator: OPFModel

function test_opf_pm(OPF::Type{<:PM.AbstractPowerModel}, casename::String)
    data = make_basic_network(pglib(casename))
    @testset "Full data" begin test_opf_pm(OPF, data) end

    if casename ∈ ["pglib_opf_case14_ieee", "pglib_opf_case30_ieee"]
        # skip n-1 tests for small cases due to infeasibility
        return
    else
        non_bridges = [e for (e, b) in OPFGenerator.bridges(data) if !b]
        drop_branch = argmin(branch -> data["branch"][branch]["rate_a"], non_bridges)
        drop_gen = argmin(gen -> gen["pmax"], values(data["gen"]))

        data_drop = deepcopy(data)
        data_drop["branch"][drop_branch]["br_status"] = 0
        data_drop["gen"]["$(drop_gen["index"])"]["gen_status"] = 0

        @testset "Outage" begin test_opf_pm(OPF, data_drop) end
    end
end

"""
    test_opf_pm(OPF, data)

Build & solve opf using OPFGenerator, and compare against PowerModels implementation.
"""
function test_opf_pm(::Type{OPF}, data::Dict) where{OPF <: PM.AbstractPowerModel}
    error("""`test_opf_pm($(OPF), data)` not implemented.
    You must implement a function with the following signature:
        function test_opf_pm(::Type{OPF}, data::Dict) where{OPF <: $(OPF)}
            # unit tests ...
            return nothing
        end
    """)

    return nothing
end

include("acp.jl")
include("dcp.jl")
include("socwr.jl")
include("ed.jl")

# other tests
include("quad_obj.jl")

const PGLIB_CASES = ["14_ieee", "30_ieee", "57_ieee", "89_pegase", "118_ieee"]

@testset "OPF" begin
    @testset "$(OPF)" for OPF in OPFGenerator.SUPPORTED_OPF_MODELS
        @testset "$(casename)" for casename in PGLIB_CASES
            test_opf_pm(OPF, "pglib_opf_case$(casename)")
        end

        @testset "QuadObj" begin test_quad_obj_warn(OPF) end
    end


    @testset _test_socwr_DualSolFormat()

    @testset _test_socwr_DualFeasibility()
end

include("utils.jl")
