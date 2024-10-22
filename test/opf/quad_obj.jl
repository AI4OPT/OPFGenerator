function test_quad_obj_warn(OPF)
    data = make_basic_network(pglib("pglib_opf_case24_ieee_rts"))

    solver = OPT_SOLVERS[OPF]

    warn_msg = "Data pglib_opf_case24_ieee_rts has quadratic cost terms; those terms are being ignored"
    opf = @test_logs (:warn, warn_msg) match_mode=:any OPFGenerator.build_opf(OPF, data, solver)

    @test isa(objective_function(opf.model), JuMP.AffExpr)

    return nothing
end
