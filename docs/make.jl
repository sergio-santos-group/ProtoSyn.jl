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
        ],
        "Peptides" => [
            "peptides/index.md"
        ]
    ]
)


