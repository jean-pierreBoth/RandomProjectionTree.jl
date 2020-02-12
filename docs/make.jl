push!(LOAD_PATH, "../src/")

DOCUMENTER_DEBUG=true

using Documenter, RandomProjectionTree


makedocs(
    format = Documenter.HTML(prettyurls = false),
    sitename = "RandomProjectionTree",
    pages = Any[
        "Introduction" => "INTRO.md",
        "RPTree.jl " => "index.md",
    ]
)

