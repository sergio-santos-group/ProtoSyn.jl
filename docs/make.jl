push!(LOAD_PATH, "../src")

using Documenter, ProtoSyn, ProtoSyn.Builder, ProtoSyn.Peptides, ProtoSyn.Sugars

makedocs(
    sitename="ProtoSyn.jl",
    pages = [
        # "Home" => "index.md",
        # "Manual" => [
        #     "manual/getting-started.md"
        # ],
        # "ProtoSyn" => [
        #     "core/index.md"
        #     "core/graph.md"
        #     "core/types.md"
        # ],
        "Builder" => "core/builder.md",
        "Peptides" => [
            "peptides/index.md"
        ],
        "Sugars" => [
            "sugars/index.md"
        ],
        # "XMLRPC" => [
        #     "xmlrpc/index.md"
        # ]
    ]
)


