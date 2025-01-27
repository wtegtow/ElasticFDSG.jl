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
        "User Guide" => [
            "General Usage" => "userguide/intro.md",
            "Velocity Models" => "userguide/velmod.md",
            "Configurations" => "userguide/config.md",
        ],
        "Method" => "method.md", 
        "API Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/wtegtow/ElasticFDSG.jl",
    devbranch="main",
)
