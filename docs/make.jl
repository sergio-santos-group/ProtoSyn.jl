push!(LOAD_PATH, "../src")

using Documenter, ProtoSyn, ProtoSyn.Peptides

DocMeta.setdocmeta!(ProtoSyn, :DocTestSetup, :(begin
    using ProtoSyn;
    using Random; Random.seed!(1);
    res_lib  = ProtoSyn.Peptides.grammar(Float64, verbose = false);
    frag     = fragment(res_lib, seq"AAA");
    pose     = build(res_lib, seq"SESEAEFKQRLAAIKTRLQAL"); sync!(pose);
    pose_mod = copy(pose);
    ProtoSyn.setdihedral!(pose_mod.state, pose_mod.graph[1][2]["N"], deg2rad(90));
    sync!(pose_mod);
    rrbm = ProtoSyn.Mutators.RotationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, randn, ProtoSyn.center_of_mass, 0.4, an"CA" & an"CB")
    trbm = ProtoSyn.Mutators.TranslationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, 1.0, nothing)
    vl = ProtoSyn.Calculators.VerletList(pose.state.size); vl.cutoff = 4.0
    ProtoSyn.Calculators.update!(ProtoSyn.SIMD_1, vl, pose)
    driver_state = ProtoSyn.Drivers.MonteCarloState{Float64}()
    dihedral_mutator = ProtoSyn.Mutators.DihedralMutator(randn, 0.01, 0.5, an"C|N"r)
    energy_function = ProtoSyn.Common.default_energy_function()
end); recursive=true)

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
                    "Monte Carlo" => "protosyn-api/core/drivers/drivers-monte-carlo.md",
                    "Steepest Descent" => "protosyn-api/core/drivers/drivers-steepest-descent.md",
                    "ILS" => "protosyn-api/core/drivers/drivers-ils.md",
                    "Compound Driver" => "protosyn-api/core/drivers/drivers-compound.md",
                ],
                "Submodules" => [
                    "Selections" => "protosyn-api/core/submodules/selections.md",
                    "Builder" => "protosyn-api/core/submodules/builder.md"
                ],
            ],
            "Peptides" => [
                "Introduction" => "protosyn-api/peptides/introduction.md",
                "Types" => "protosyn-api/peptides/types.md"
            ]
        ]
    ],
    doctest = true,
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
