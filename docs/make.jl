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
            "I/O utilities" => "io.md",
        ],
        "Reference" => "lib/public.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(repo="github.com/AI4OPT/OPFGenerator.git")
