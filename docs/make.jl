using OPTinv
#using UnitfulLinearAlgebra # incl b.c. non-registered
using TMI 
using Documenter

DocMeta.setdocmeta!(OPTinv, :DocTestSetup, :(using OPTinv); recursive=true)

makedocs(;
    modules=[OPTinv],
    authors="B. Hamilton <brynnydd.hamilton@whoi.edu>",
    repo="https://github.com/b-r-hamilton/OPTinv.jl/blob/{commit}{path}#{line}",
    sitename="OPTinv.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://b-r-hamilton.github.io/OPTinv.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/b-r-hamilton/OPTinv.jl",
    devbranch="main",
)
