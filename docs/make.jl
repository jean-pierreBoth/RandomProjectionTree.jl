
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

deploydocs(
    repo = "RandomProjectionTree.git",
    target = "build",
    julia = "1.3",
    osname = "linux",
    deps = nothing,
    make = nothing
)
