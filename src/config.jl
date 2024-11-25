"""
    _get_case_info(config)

Extract case file and name from input `config`.

To be valid, the input config should include:
* either a `case_file` or `pglib_case` entry
* if no `case_file` is provided, `pglib_case` should be a valid, unique PGLib case name.

The case name will be set to the generic \"case\" value if none is provided.

!!! warning
    if `case_file` is provided, `pglib_case` will be ignored.
    Therefore, users should provide `case_name` when supplying `case_file`.
"""
function _get_case_info(config)
    case_file = if haskey(config, "case_file")
        config["case_file"]
    elseif haskey(config, "pglib_case")
        # Check that case reference is valid PGLib name
        pglib_cases = PGLib.find_pglib_case(config["pglib_case"])
        if length(pglib_cases) == 1
            # All good
            joinpath(PGLib.PGLib_opf, pglib_cases[1])
        elseif length(pglib_cases) == 0
            error("PGLib case `$(case_name)` was not found. Try running `PGLib.find_pglib_case(\"$(case_name)\")` to find similar case names.")
        else
            error("""PGLib case returned matches; please update the case name to be more specific.
            Matching PGLib cases:\n""" * prod("* $(case)\n" for case in pglib_cases))
        end
    else
        error("Invalid config: must provide either \"case_file\" or \"pglib_case\".")
    end

    case_name = splitext(basename(case_file))[1]
    case_name == "" && (case_name = "case")  # generic fallback, just in case

    return case_file, case_name
end
