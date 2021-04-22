push!(LOAD_PATH, "../src")

using Documenter, ProtoSyn, ProtoSyn.Peptides

makedocs(
    sitename="ProtoSyn.jl",
    authors="José Pereira & Sérgio Santos",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
    ),
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
                "Methods" => [
                    "protosyn-api/core/methods/graph.md",
                    "protosyn-api/core/methods/state.md",
                    "protosyn-api/core/methods/pose.md",
                    "protosyn-api/core/methods/io.md",
                    "protosyn-api/core/methods/aux.md",
                ],
                "Submodules" => [
                    "Selections" => "protosyn-api/core/submodules/selections.md",
                    "Builder" => "protosyn-api/core/submodules/builder.md"
                ]
            ],
            "Peptides" => [
                "Introduction" => "protosyn-api/peptides/introduction.md"
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
