using Documenter
using OPFGenerator

makedocs(
    modules=[OPFGenerator],
    sitename = "OPFGenerator",
    format = Documenter.HTML(;
        mathengine = Documenter.MathJax(),
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Problem formulations" => [
                "Notations" => "opf/notations.md",
                "AC-OPF"    => "opf/acp.md",
                "SOC-OPF"   => "opf/socwr.md",
                "DC-OPF"    => "opf/dcp.md",
            ],
            "I/O utilities" => "io.md",
        ],
        "Reference" => "lib/public.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo="github.com/AI4OPT/OPFGenerator.git",
    push_preview=true,
)
