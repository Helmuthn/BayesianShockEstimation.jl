using BayesianShockEstimation
using Documenter

DocMeta.setdocmeta!(BayesianShockEstimation, :DocTestSetup, :(using BayesianShockEstimation); recursive=true)

makedocs(;
    modules=[BayesianShockEstimation],
    authors="Helmuth Naumer <hnaumer2@illinois.edu> and contributors",
    repo="https://github.com/Helmuthn/BayesianShockEstimation.jl/blob/{commit}{path}#{line}",
    sitename="BayesianShockEstimation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Helmuthn.github.io/BayesianShockEstimation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Helmuthn/BayesianShockEstimation.jl",
    devbranch="main",
)
