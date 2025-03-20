var documenterSearchIndex = {"docs":
[{"location":"opf/notations/#Notations","page":"Notations","title":"Notations","text":"","category":"section"},{"location":"opf/notations/#Input-parameters","page":"Notations","title":"Input parameters","text":"","category":"section"},{"location":"opf/notations/","page":"Notations","title":"Notations","text":"Sets:","category":"page"},{"location":"opf/notations/","page":"Notations","title":"Notations","text":"mathcalN = 1  N: Set of buses\nmathcalE = 1  E: Set of branches\nmathcalE^+_i: Set of branches leaving bus i in mathcalN\nmathcalE^-_i: Set of branches entering bus i in mathcalN\nmathcalG = 1  G: Set of generators\nmathcalG_i = 1  G_i: Set of generators at bus i in mathcalN\nmathcalL = 1  L: Set of loads at bus i in mathcalN\nmathcalS = 1  S: Set of shunts\nmathcalS = 1  S_i: Set of shunts at bus i in mathcalN.","category":"page"},{"location":"opf/notations/","page":"Notations","title":"Notations","text":"Network data:","category":"page"},{"location":"opf/notations/","page":"Notations","title":"Notations","text":"g^s_s + mathbfj  b^s_s: complex admittance of shunt s in mathcalS\ng_e + mathbfj  b_e: complex admittance of branch e in mathcalE","category":"page"},{"location":"opf/dcp/#DC-OPF","page":"DC-OPF","title":"DC-OPF","text":"","category":"section"},{"location":"opf/dcp/","page":"DC-OPF","title":"DC-OPF","text":"See DCPPowerModel formulation in PowerModels.jl.","category":"page"},{"location":"opf/dcp/#Mathematical-formulation","page":"DC-OPF","title":"Mathematical formulation","text":"","category":"section"},{"location":"opf/dcp/","page":"DC-OPF","title":"DC-OPF","text":"The primal problem reads","category":"page"},{"location":"opf/dcp/","page":"DC-OPF","title":"DC-OPF","text":"beginalign\n    min_mathbfp^textg mathbfpf boldsymboltheta quad \n     labeleqDCPobjective\n        sum_g in mathcalG c_g mathbfp^textg_g + c^0_g\n        \n    textst quad\n     labeleqDCPslack_bus\n        boldsymboltheta_textslack = 0\n        \n         lambda^textslack\n     labeleqDCPkirchhoff\n        sum_g in mathcalG_i mathbfp^textg_g \n        - sum_e in mathcalE^+_i mathbfpf_e\n        + sum_e in mathcalE^-_i mathbfpf_e\n        = \n        sum_l in mathcalL_i p^d_l\n        + sum_s in mathcalS_i g^s_s\n         forall i in mathcalN\n         lambda^textkcl\n     labeleqDCPohm\n        mathbfpf_e\n        =\n        b_e (boldsymboltheta_t(e) - boldsymboltheta_s(e))\n         forall e in mathcalE\n         lambda^textohm\n     labeleqDCPva_diff\n        underlineDelta boldsymboltheta_e\n        leq\n        boldsymboltheta_t(e) - boldsymboltheta_s(e)\n        leq \n        overlineDelta boldsymboltheta_e\n         forall e in mathcalE\n         mu^Delta boldsymboltheta\n     labeleqDCPthermal\n        -overlinef_e leq mathbfpf_e leq overlinef_e\n         forall e in mathcalE\n         mu^f\n     labeleqDCPpg_bounds\n        underlinep^g_g leq mathbfp^textg_g leq overlinep^g_g\n         forall g in mathcalG\n         mu^pg\nendalign","category":"page"},{"location":"opf/dcp/#Nomenclature","page":"DC-OPF","title":"Nomenclature","text":"","category":"section"},{"location":"opf/dcp/#Primal-variables","page":"DC-OPF","title":"Primal variables","text":"","category":"section"},{"location":"opf/dcp/","page":"DC-OPF","title":"DC-OPF","text":"Symbol Data Size Description\nmathbfp^textg pg G Active power dispatch\nboldsymboltheta va N Nodal voltage angle\nmathbfpf pf E Active power flow","category":"page"},{"location":"opf/dcp/#Dual-variables","page":"DC-OPF","title":"Dual variables","text":"","category":"section"},{"location":"opf/dcp/","page":"DC-OPF","title":"DC-OPF","text":"Symbol Data Size Associated constraint\nlambda^textkcl lam_kirchhoff N Nodal power balance eqrefeqDCPkirchhoff\nlambda^textohm lam_ohm E Ohm's law eqrefeqDCPohm\nmu^Delta boldsymboltheta mu_va_diff E Angle difference limit eqrefeqDCPva_diff\nmu^pf mu_sm E Thermal limit eqrefeqDCPthermal\nmu^pg mu_pg G Generation min/max limits eqrefeqDCPpg_bounds\nlambda^textslack – 1 Slack bus voltage angle eqrefeqDCPslack_bus","category":"page"},{"location":"opf/dcp/","page":"DC-OPF","title":"DC-OPF","text":"Dual variable lambda^textslack is always zero at the optimum, hence it is not exported.","category":"page"},{"location":"opf/acp/#AC-OPF","page":"AC-OPF","title":"AC-OPF","text":"","category":"section"},{"location":"opf/acp/","page":"AC-OPF","title":"AC-OPF","text":"See ACPPowerModel formulation in PowerModels.jl.","category":"page"},{"location":"lib/public/","page":"Reference","title":"Reference","text":"Modules = [OPFGenerator]","category":"page"},{"location":"lib/public/#OPFGenerator.E2ELRReserveScaler","page":"Reference","title":"OPFGenerator.E2ELRReserveScaler","text":"E2ELRReserveScaler\n\nSamples reserve requirements following the procedure below:\n\nSample a minimum reserve requirement MRR from a uniform distribution U(lb, ub) (mrr_dist).\nCompute the upper bound of reserve requirements for each generator as rmax = α * (pmax - pmin).\nFix the lower bound of reserve requirement per generator to zero.\nFix the reserve cost of each generator to zero.\n\nThe parameter α is a scaling factor that determines each generator's maximum reserve     capacity. It is the factor parameter times the ratio of the largest generator's capacity     to the sum of all generators' dispatchable capacity.\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#OPFGenerator.Glocal","page":"Reference","title":"OPFGenerator.Glocal","text":"Glocal{G,L}\n\nA glocal distribution with global/local factors α::G and  η::L.\n\nThis distribution represents a random variable of the form ϵ = α×η, where\n\nα is a scalar random variable, with distribution d_α::G\nη is a vector random variable, with distribution d_η::L\nα and η are independent random variables\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#OPFGenerator.LoadScaler","page":"Reference","title":"OPFGenerator.LoadScaler","text":"LoadScaler{D}\n\nScales loads with multiplicative noise sample from d::D.\n\nScaling retains power factors, i.e., each load's active      and reactive demand is scaled by the same number.\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#OPFGenerator.ScaledLogNormal-Tuple{Float64, Float64, Vector{Float64}}","page":"Reference","title":"OPFGenerator.ScaledLogNormal","text":"ScaledLogNormal(l, u, σs)\n\nGenerate a Glocal distribution ϵ = α×η where α~U[l,u] and ηᵢ ~ LogNormal(-σᵢ²/2, σᵢ).\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.ScaledUniform-Tuple{Float64, Float64, Vector{Float64}}","page":"Reference","title":"OPFGenerator.ScaledUniform","text":"ScaledUniform(l, u, σs)\n\nGenerate a Glocal distribution ϵ = α×η where α ~ U[l,u] and ηᵢ ~ U[1-σᵢ, 1+σᵢ].\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator._dedupe_and_sort_h5!-Tuple{Any}","page":"Reference","title":"OPFGenerator._dedupe_and_sort_h5!","text":"_dedupe_and_sort_h5!(D)\n\nDe-duplicated and sort dataset D in increasing order of random seeds.\n\nEquivalent to _dedupe_h5!(D); _sort_h5!(D), but more efficient.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator._dedupe_h5!-Tuple{Any}","page":"Reference","title":"OPFGenerator._dedupe_h5!","text":"_dedupe_h5!(D)\n\nDe-duplicate points in h5 dataset D, according to their random seed.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator._merge_h5-Union{Tuple{Array{Array{T, N}, 1}}, Tuple{N}, Tuple{T}} where {T, N}","page":"Reference","title":"OPFGenerator._merge_h5","text":"_merge_h5(V::Vector{Array{T,N})\n\nConcatenate a collection of N-dimensional arrays along their last dimension.\n\nThis function is semantically equivalent to cat(V...; dims=ndims(first(V))),     but uses a more efficient, splatting-free, implementation. All elements of V must have the same size in the first N-1 dimensions.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator._select_h5!-Tuple{Dict, Any}","page":"Reference","title":"OPFGenerator._select_h5!","text":"_select_h5!(D, p)\n\nSelect data points in D as indicated by p.\n\nD should be a dictionary in h5-compatible format, and p is either a     vector of indices, or a logical vector of the same length as D[\"meta\"][\"seed\"].\n\nIf p is a vector of indices, then all values of p should be integers    between 1 and the number of elements in D\nIf p is a logical vector, then it should have the same length as D[\"meta\"][\"seed\"].   Only datapoints i for which p[i] is true are selected.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator._sort_h5!-Tuple{Any}","page":"Reference","title":"OPFGenerator._sort_h5!","text":"_sort_h5!(D)\n\nSort dataset D in increasing order of random seeds.\n\nThe dictionary D should be in h5-compatible format. It is modified in-place.\n\nThe function expects D[\"meta\"][\"seed\"] to exist and be a Vector{Int}.     An error is thrown if such an entry is not found.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.bridges-Tuple{Dict}","page":"Reference","title":"OPFGenerator.bridges","text":"bridges(data)\n\nIdentify whether each branch is a bridge.\n\nThe input data must be in basic format.\n\nA branch is a bridge if removing it renders the network disconnected. Returns a dictionary res::Dict{String,Bool} such that     res[br] is true if branch br is a bridge, and false otherwise.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.build_opf-Tuple{Type{PowerModels.ACPPowerModel}, Dict{String, Any}, Any}","page":"Reference","title":"OPFGenerator.build_opf","text":"build_acopf(data, optimizer)\n\nBuild an AC-OPF model.\n\nThis implementation is based on the AC-OPF formulation of Rosetta-OPF     https://github.com/lanl-ansi/rosetta-opf/blob/38a951326df3156d79dcdc49c8010aa29905b05d/jump.jl\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.build_opf-Tuple{Type{PowerModels.DCPPowerModel}, Dict{String, Any}, Any}","page":"Reference","title":"OPFGenerator.build_opf","text":"build_dcopf(data, optimizer)\n\nBuild a DC-OPF model.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.build_opf-Union{Tuple{OPF}, Tuple{Type{OPF}, Dict{String, Any}, Any}} where OPF<:Union{PowerModels.SOCWRConicPowerModel, PowerModels.SOCWRPowerModel}","page":"Reference","title":"OPFGenerator.build_opf","text":"build_soc_opf(data, optimizer)\n\nBuild an SOC-OPF model.\n\nThis implementation is based on the SOC-OPF formulation of PMAnnex.jl     https://github.com/lanl-ansi/PMAnnex.jl/blob/f303f3c3c61e2d1a050ee7651fa6e8abc4055b55/src/model/opf.jl\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.extract_result-Tuple{OPFGenerator.OPFModel{PowerModels.ACPPowerModel}}","page":"Reference","title":"OPFGenerator.extract_result","text":"_extract_acopf_solution(model, data)\n\nExtract ACOPF solution from optimization model. The model must have been solved before.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.extract_result-Tuple{OPFGenerator.OPFModel{PowerModels.DCPPowerModel}}","page":"Reference","title":"OPFGenerator.extract_result","text":"_extract_dcopf_solution(model, data)\n\nExtract DCOPF solution from optimization model. The model must have been solved before.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.extract_result-Union{Tuple{OPFGenerator.OPFModel{OPF}}, Tuple{OPF}} where OPF<:Union{PowerModels.SOCWRConicPowerModel, PowerModels.SOCWRPowerModel}","page":"Reference","title":"OPFGenerator.extract_result","text":"_extract_solution(model, data)\n\nExtract SOC-OPF solution from optimization model. The model must have been solved before.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.load_h5","page":"Reference","title":"OPFGenerator.load_h5","text":"load_h5\n\n\n\n\n\n","category":"function"},{"location":"lib/public/#OPFGenerator.load_json-Tuple{AbstractString}","page":"Reference","title":"OPFGenerator.load_json","text":"load_json(filename::AbstractString)\n\nLoad JSON data from file filename.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.save_h5-Tuple{AbstractString, Any}","page":"Reference","title":"OPFGenerator.save_h5","text":"save_h5(filename, D)\n\nSaves dictionary D to HDF5 file filename.\n\nAll keys in D must be of String type, and it must be HDF5-compatible. Additional restrictions are enforced on the values of D, see below.\n\nwarning: Warning\nOnly the following types are supported:String\n(un)signed integers up to 64-bit precision\nFloat32 and Float64\nComplex versions of the above numeric types\nDense Arrays of the the above scalar typesNumerical data of an unsupported type will be converted to Float64 when possible. An error will be thrown if the conversion is not possible.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.save_json-Tuple{AbstractString, Any}","page":"Reference","title":"OPFGenerator.save_json","text":"save_json(filename::AbstractString, data; indent)\n\nSave data into JSON file filename. The following formats are supported:\n\nuncompressed JSON .json\nGzip-compressed JSON .json.gz\nBzip2-compressed JSON .json.bz2\n\nIf the file extension does not match one of the above, an error is thrown.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#OPFGenerator.tensorize-Union{Tuple{Vector{T}}, Tuple{T}} where T","page":"Reference","title":"OPFGenerator.tensorize","text":"tensorize(V)\n\nConcatenate elements of V into a higher-dimensional tensor.\n\nSimilar to Base.stack, with one major difference: if V is a vector of scalars,     the result is a 2D array M whose last dimension is length(V),     and such that M[:, i] == V[i].\n\nThis function is only defined for Vector{T} and Vector{Array{T,N}} inputs,     to avoid any unexpected behavior of Base.stack.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#Random.rand!-Tuple{Random.AbstractRNG, SimpleOPFSampler, Dict}","page":"Reference","title":"Random.rand!","text":"rand!(rng::AbstractRNG, s::AbstractOPFSampler, data::Dict)\n\nSample one new OPF instance and modify data in-place.\n\ndata must be a Dict in PowerModels format, representing the same network     (i.e., same grid components with same indexing) as the one used to create s.\n\n\n\n\n\n","category":"method"},{"location":"#OPFGenerator.jl","page":"Home","title":"OPFGenerator.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for OPFGenerator.jl","category":"page"},{"location":"","page":"Home","title":"Home","text":"warning: Warning\nThis documentation is a work in progress. Please open an issue if content is missing / erroneous","category":"page"},{"location":"io/","page":"I/O utilities","title":"I/O utilities","text":"OPFGenerator provides utilities to save/load files in the HDF5 and JSON formats.","category":"page"},{"location":"io/#HDF5","page":"I/O utilities","title":"HDF5","text":"","category":"section"},{"location":"io/","page":"I/O utilities","title":"I/O utilities","text":"To load HDF5 files, use the load_h5, which is the same as HDF5.h5read.","category":"page"},{"location":"io/","page":"I/O utilities","title":"I/O utilities","text":"h = load_h5(\"test_file.h5\", \"/\")  # load all HDF5 file","category":"page"},{"location":"io/","page":"I/O utilities","title":"I/O utilities","text":"To save HDF5 files to disk, use the save_h5 function. All keys of the dictionary (and sub-dictionaries) must be Strings.","category":"page"},{"location":"io/","page":"I/O utilities","title":"I/O utilities","text":"save_h5(\"myfile.h5\", d)","category":"page"},{"location":"io/","page":"I/O utilities","title":"I/O utilities","text":"warning: Warning\nWhen exporting data to an HDF5 file, only the following data types are supported:String\nNumber types supported by HDF5.jl\nArrays of thoseNumerical data in an unsupported type will be converted to Float64, when possible.","category":"page"},{"location":"io/#JSON","page":"I/O utilities","title":"JSON","text":"","category":"section"},{"location":"io/","page":"I/O utilities","title":"I/O utilities","text":"warning: Warning\nOnly .json extensions are supported","category":"page"},{"location":"io/","page":"I/O utilities","title":"I/O utilities","text":"Use the load_json and save_json functions to load/save data to/from JSON files.","category":"page"},{"location":"io/","page":"I/O utilities","title":"I/O utilities","text":"using OPFGenerator\n\n# Load a dictionary from a JSON file\nd = load_json(\"my_json_file.json\")\n\n# Save a dictionary to a JSON file\nsave_json(\"my_new_jsonfile.json\", d)\nsave_json(\"my_pretty_jsonfile.json\", d, indent=2)  # prettier formatting","category":"page"},{"location":"io/","page":"I/O utilities","title":"I/O utilities","text":"Compressed JSON files (.json.gz and .json.bz2) are supported automatically","category":"page"},{"location":"io/","page":"I/O utilities","title":"I/O utilities","text":"# Load a dictionary from a compressed JSON file\nd = load_json(\"my_json_file.json.gz\")\nd = load_json(\"my_json_file.json.bz2\")\n\n# Save a dictionary to a compressed JSON file\nsave_json(\"my_new_jsonfile.json.gz\", d)\nsave_json(\"my_new_jsonfile.json.bz2\", d)","category":"page"},{"location":"opf/socwr/#SOC-OPF","page":"SOC-OPF","title":"SOC-OPF","text":"","category":"section"},{"location":"opf/socwr/","page":"SOC-OPF","title":"SOC-OPF","text":"See SOCWRPowerModel and SOCWRConicPowerModel formulations in PowerModels.jl.","category":"page"},{"location":"opf/socwr/","page":"SOC-OPF","title":"SOC-OPF","text":"warning: Warning\nThe formulation in OPFGenerator is not equivalent to that implemented in PowerModels","category":"page"}]
}
