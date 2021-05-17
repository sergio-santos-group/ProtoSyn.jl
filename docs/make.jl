push!(LOAD_PATH, "../src")

using Documenter, ProtoSyn, ProtoSyn.Peptides

DocMeta.setdocmeta!(ProtoSyn, :DocTestSetup, :(using ProtoSyn); recursive=true)

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
                "Calculators" => [
                    "Calculators Section" => "protosyn-api/core/calculators/calculators-section.md",
                    "TorchANI" => "protosyn-api/core/calculators/torchani.md",
                    "Bond distance Restraint" => "protosyn-api/core/calculators/bond-distance-restraint.md",
                    "Potential Restraints" => "protosyn-api/core/calculators/potential-restraints.md",
                ],
                "Mutators" => [
                    "Mutators Section" => "protosyn-api/core/mutators/mutators-section.md",
                    "Dihedral Mutator" => "protosyn-api/core/mutators/mutators-dihedral.md",
                    "Crankshaft Mutator" => "protosyn-api/core/mutators/mutators-crankshaft.md",
                    "Rigid Body Mutators" => "protosyn-api/core/mutators/mutators-rigid-body.md",
                    "Compound Mutators" => "protosyn-api/core/mutators/mutators-compound.md",
                ],
                "Drivers" => [
                    "Drivers Section" => "protosyn-api/core/drivers/drivers-section.md",
                ],
                "Submodules" => [
                    "Selections" => "protosyn-api/core/submodules/selections.md",
                    "Builder" => "protosyn-api/core/submodules/builder.md"
                ],
            ],
            "Peptides" => [
                "Introduction" => "protosyn-api/peptides/introduction.md"
            ]
        ]
    ],
    doctest = false,
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
