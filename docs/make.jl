using PopSim
using Documenter

DocMeta.setdocmeta!(PopSim, :DocTestSetup, :(using PopSim); recursive=true)

makedocs(;
    modules=[PopSim],
    authors="Peter Arndt <arndt@molgen.mpg.de> and contributors",
    sitename="PopSim.jl",
    format=Documenter.HTML(;
        canonical="https://ArndtLab.github.io/PopSim.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    warnonly=[:missing_docs],
)

deploydocs(;
    repo="github.com/ArndtLab/PopSim.jl",
    devbranch="main",
)
