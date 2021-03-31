push!(LOAD_PATH, "../src")

using Documenter, ProtoSyn, ProtoSyn.Builder, ProtoSyn.Peptides

makedocs(
    sitename="ProtoSyn.jl",
    modules=[ProtoSyn],
    pages = [
        "Home" => "index.md",
        "Getting Started" => [
            "Installation" => "getting-started/installation.md",
            "First steps" => "getting-started/first-steps.md"
        ],
        "ProtoSyn API" => [
            "Core" => [
                "Types" => "protosyn-api/core/types.md",
                "Selections" => "protosyn-api/core/selections.md"
            ]
        ]
    ]
)

deploydocs(
    repo = "github.com/sergio-santos-group/ProtoSyn.jl.git",
    # osname = "linux",
    # julia = "1.5",
    #deps = nothing,
    #make = nothing,
    #target = "build",
    branch = "use-cases",
    devurl = "dev",
)
