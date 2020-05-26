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

deploydocs(
    repo = "github.com/sergio-santos-group/ProtoSyn.jl.git",
    # osname = "linux",
    # julia = "1.0",
    #deps = nothing,
    #make = nothing,
    #target = "build",
    #devbranch = "refactoring",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#", "dev" => "dev", "v0.2-alpha"=>"v0.2-alpha"],
)
