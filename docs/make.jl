using ElasticFDSG, ElasticFDSG.dim2, ElasticFDSG.dim3
using Documenter

DocMeta.setdocmeta!(ElasticFDSG, :DocTestSetup, :(using ElasticFDSG, ElasticFDSG.dim2, ElasticFDSG.dim3); recursive=true)

makedocs(;
    modules=[ElasticFDSG, ElasticFDSG.dim2, ElasticFDSG.dim3],
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
        "Validation" => "validate.md",
        "API Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/wtegtow/ElasticFDSG.jl",
    devbranch="main",
)
