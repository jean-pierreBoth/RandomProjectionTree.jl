push!(LOAD_PATH, "../src/")

DOCUMENTER_DEBUG=true

using Documenter, RandomProjectionTree


makedocs(
    format = :html,
    sitename = "RandomProjectionTree",
    pages = Any[
        "Introduction" => "INTRO.md",
        "RPTree.jl documentation" => "index.md",
    ]
)

