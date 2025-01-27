using ElasticFDSG
using Documenter

DocMeta.setdocmeta!(ElasticFDSG, :DocTestSetup, :(using ElasticFDSG); recursive=true)

makedocs(;
    modules=[ElasticFDSG],
    authors="William Tegtow <w.tegtow@gmail.com> and contributors",
    sitename="ElasticFDSG.jl",
    format=Documenter.HTML(;
        canonical="https://wtegtow.github.io/ElasticFDSG.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/wtegtow/ElasticFDSG.jl",
    devbranch="main",
)
