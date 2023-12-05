using VaidyaPT
using Documenter

DocMeta.setdocmeta!(VaidyaPT, :DocTestSetup, :(using VaidyaPT); recursive=true)

makedocs(;
    modules=[VaidyaPT],
    authors="Jaime <jaime.redondo.yuste@gmail.com> and contributors",
    repo="https://github.com/jredondoyuste/VaidyaPT.jl/blob/{commit}{path}#{line}",
    sitename="VaidyaPT.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jredondoyuste.github.io/VaidyaPT.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jredondoyuste/VaidyaPT.jl",
    devbranch="main",
)
