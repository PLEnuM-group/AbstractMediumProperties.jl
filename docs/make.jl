using AbstractMediumProperties
using Documenter

DocMeta.setdocmeta!(AbstractMediumProperties, :DocTestSetup, :(using AbstractMediumProperties); recursive=true)

makedocs(;
    modules=[AbstractMediumProperties],
    authors="Christian Haack <chr.hck@gmail.com>",
    sitename="AbstractMediumProperties.jl",
    format=Documenter.HTML(;
        canonical="https://chrhck.github.io/AbstractMediumProperties.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "Guide" => "guide.md",
    ],
    nav=[
        "Home" => "index.md",
        "API" => "api.md",
        "Guide" => "guide.md",
    ],
)

deploydocs(;
    repo="github.com/PLEnuM-group/AbstractMediumProperties.jl",
    devbranch="main",
)
