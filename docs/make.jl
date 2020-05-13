push!(LOAD_PATH, "../src")

using Documenter, ProtoSyn, ProtoSyn.Peptides

makedocs(
    sitename="ProtoSyn.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "manual/getting-started.md"
        ],
        "Core" => [
            "core/index.md"
            "core/graph.md"
        ],
        "Peptides" => [
            "peptides/index.md"
        ]
    ]
)


