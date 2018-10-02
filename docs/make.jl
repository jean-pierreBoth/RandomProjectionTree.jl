push!(LOAD_PATH, "/media/jpboth/scandisk1/Julia/")

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
    julia = "0.6",
    osname = "linux",
    deps = nothing,
    make = nothing
)
