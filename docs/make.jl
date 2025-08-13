using APop
using Documenter

DocMeta.setdocmeta!(APop, :DocTestSetup, :(using APop); recursive=true)

makedocs(;
    modules=[APop],
    authors="Peter Arndt <arndt@molgen.mpg.de> and contributors",
    sitename="APop.jl",
    format=Documenter.HTML(;
        canonical="https://ArndtLab.github.io/APop.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    warnonly=[:missing_docs],
)

deploydocs(;
    repo="github.com/ArndtLab/APop.jl",
    devbranch="main",
)
