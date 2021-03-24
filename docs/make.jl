using OctreeBH
using Documenter

DocMeta.setdocmeta!(OctreeBH, :DocTestSetup, :(using OctreeBH); recursive=true)

makedocs(;
    modules=[OctreeBH],
    authors="Chia-Yu Hu <cyhu.astro@gmail.com> and contributors",
    repo="https://github.com/huchiayu/OctreeBH.jl/blob/{commit}{path}#{line}",
    sitename="OctreeBH.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://huchiayu.github.io/OctreeBH.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/huchiayu/OctreeBH.jl",
)
