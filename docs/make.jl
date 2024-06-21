using Documenter
using OPFGenerator

makedocs(
    sitename = "OPFGenerator",
    format = Documenter.HTML(;
        mathengine = Documenter.MathJax(),
    ),
    pages = [
        "Home" => "index.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(repo="github.com/AI4OPT/OPFGenerator.git")
