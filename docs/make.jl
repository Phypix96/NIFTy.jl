using NIFTy
using Documenter

makedocs(;
    modules=[NIFTy],
    authors="Philipp Haim <philipp.haim@gmail.com> and contributors",
    repo="https://github.com/Philipp Haim/NIFTy.jl/blob/{commit}{path}#L{line}",
    sitename="NIFTy.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Philipp Haim.gitlab.io/NIFTy.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
