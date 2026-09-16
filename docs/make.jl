using GraphsHDC
using Documenter

DocMeta.setdocmeta!(GraphsHDC, :DocTestSetup, :(using GraphsHDC); recursive = true)

makedocs(;
    modules = [GraphsHDC],
    authors = "Carlos Vigil-Vásquez",
    repo = "https://github.com/cvigilv/GraphsHDC.jl/blob/{commit}{path}#{line}",
    sitename = "GraphsHDC.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://cvigilv.github.io/KEGGAPI.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "man/getting-started.md",
        # "Examples" => Any[ ],
        "API reference" => "man/api.md",
    ],
)

deploydocs(;
    repo = "github.com/cvigilv/GraphsHDC.jl",
    devbranch = "main",
)

