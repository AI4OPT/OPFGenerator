using OPFGenerator: OPFModel

function test_opf_pm(OPF::Type{<:PM.AbstractPowerModel}, casename::String)
    data = make_basic_network(pglib(casename))
    return test_opf_pm(OPF, data)
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
include("sdpwrm.jl")

# other tests
include("quad_obj.jl")

const PGLIB_CASES = ["14_ieee", "30_ieee", "57_ieee", "89_pegase", "118_ieee"]
const PGLIB_SDP_CASES = ["5_pjm", "14_ieee"]

@testset "OPF" begin
    @testset "$(OPF)" for OPF in OPFGenerator.SUPPORTED_OPF_MODELS
        if OPF == PowerModels.SDPWRMPowerModel
            cases = PGLIB_SDP_CASES
        else
            cases = PGLIB_CASES
        end
        @testset "$(casename)" for casename in cases
            test_opf_pm(OPF, "pglib_opf_case$(casename)")
        end

        @testset "QuadObj" begin test_quad_obj_warn(OPF) end
    end


    @testset _test_socwr_DualSolFormat()

    @testset _test_socwr_DualFeasibility()
end
