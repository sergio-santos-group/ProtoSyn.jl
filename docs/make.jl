push!(LOAD_PATH, "../src")

using Documenter, ProtoSyn, ProtoSyn.Peptides
ProtoSyn.set_logger_to_warn()

DocMeta.setdocmeta!(ProtoSyn, :DocTestSetup, :(begin
    using ProtoSyn;
    using ProtoSyn.Units;
    using Random; Random.seed!(1);
    res_lib  = ProtoSyn.Peptides.grammar;
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
    monte_carlo = ProtoSyn.Drivers.MonteCarlo(energy_function, dihedral_mutator, nothing, 10, ProtoSyn.Drivers.get_linear_quench(1.0, 10))
    cb = ProtoSyn.Common.default_energy_step_frame_callback(10, "teste.pdb")
end); recursive=true)


makedocs(
    sitename="ProtoSyn.jl",
    authors="José Pereira & Sérgio Santos",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        ),
    doctest = false,
    checkdocs=:exports,
    modules=[ProtoSyn],
    pages = [
        "Home" => "index.md",
        "Getting Started" => [
            "Installation" => "getting-started/installation.md",
            "First steps" => "getting-started/first-steps.md",
            "Examples" => "getting-started/examples.md"
        ],
        "ProtoSyn API" => [
            "Core" => [
                "Types" => "protosyn-api/core/types.md",
                "Methods" => [
                    "protosyn-api/core/methods/graph.md",
                    "protosyn-api/core/methods/state.md",
                    "protosyn-api/core/methods/pose.md",
                    "protosyn-api/core/methods/io.md",
                    "protosyn-api/core/methods/other.md",
                ],
                "Calculators" => [
                    "Calculators Section" => "protosyn-api/core/calculators/calculators-section.md",
                    "TorchANI" => "protosyn-api/core/calculators/torchani.md",
                    "Bond distance Restraint" => "protosyn-api/core/calculators/bond-distance-restraint.md",
                    "Potential Restraints" => "protosyn-api/core/calculators/potential-restraints.md",
                    "Electrostatics" => "protosyn-api/core/calculators/electrostatics.md",
                    "Generalized Born" => "protosyn-api/core/calculators/gb-solvation.md",
                    "SASA" => "protosyn-api/core/calculators/sasa.md",
                    "Hydrogen Bonds" => "protosyn-api/core/calculators/hydrogen-bonds.md",
                    "Radius of gyration" => "protosyn-api/core/calculators/radius-gyration.md",
                    "Custom reference energy" => "protosyn-api/core/calculators/custom-ref-energy.md",
                    "REF-15" => "protosyn-api/core/calculators/ref15.md",
                ],
                "Mutators" => [
                    "Mutators Section" => "protosyn-api/core/mutators/mutators-section.md",
                    "Dihedral Mutator" => "protosyn-api/core/mutators/mutators-dihedral.md",
                    "Crankshaft Mutator" => "protosyn-api/core/mutators/mutators-crankshaft.md",
                    "Rigid Body Mutators" => "protosyn-api/core/mutators/mutators-rigid-body.md",
                    "Backrub Mutators" => "protosyn-api/core/mutators/mutators-backrub.md",
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
                    "Builder" => "protosyn-api/core/submodules/builder.md",
                    "External packages" => "protosyn-api/core/submodules/external-packages.md",
                ],
            ],
            "Peptides" => [
                "Introduction" => "protosyn-api/peptides/introduction.md",
                "Types" => "protosyn-api/peptides/types.md",
                "Methods" => [
                    "Input and Output (IO)" => "protosyn-api/peptides/methods/io.md",
                    "Graph" => "protosyn-api/peptides/methods/graph.md",
                    "State" => "protosyn-api/peptides/methods/state.md",
                    "Pose" => "protosyn-api/peptides/methods/pose.md"
                ],
                "Calculators" => [
                    "Potential Restraints" => "protosyn-api/peptides/calculators/potential-restraints.md",
                    "Electrostatics" => "protosyn-api/peptides/calculators/electrostatics.md",
                    "SASA" => "protosyn-api/peptides/calculators/sasa.md",
                    "Natural frequency" => "protosyn-api/peptides/calculators/natural-frequency.md",
                    "Secondary structure propensity" => "protosyn-api/peptides/calculators/ss-propensity.md",
                    "Caterpilar solvation model" => "protosyn-api/peptides/calculators/caterpillar-solvation.md",
                    "SeqDes model" => "protosyn-api/peptides/calculators/seqdes.md",
                ],
                "Mutators" => [
                    "Rotamer Mutator" => "protosyn-api/peptides/mutators/rotamer.md",
                    "Design Mutator" => "protosyn-api/peptides/mutators/design.md"
                ],
                "Drivers" => [
                    "Rotamer Blitz" => "protosyn-api/peptides/drivers/rotamer-blitz-driver.md"
                ],
                "Submodules" => [
                    "Selections" => "protosyn-api/peptides/submodules/selections.md",
                    "Builder" => "protosyn-api/peptides/submodules/builder.md",
                    "Rotamers" => "protosyn-api/peptides/submodules/rotamers.md",
                    "External packages" => "protosyn-api/peptides/submodules/external-packages.md",
                ]
            ],
            "Materials" => [
                "Introduction" => "protosyn-api/materials/introduction.md",
                "Methods" => [
                    "Lattices" => "protosyn-api/materials/methods/lattices.md",
                    "Carbons" => "protosyn-api/materials/methods/carbons.md",
                ]
            ],
            "Sugars" => [
                "Introduction" => "protosyn-api/sugars/introduction.md",
                "Submodules" => [
                    "Builder" => "protosyn-api/sugars/submodules/builder.md",
                ]
            ],
            "Common" => [
                "Introduction" => "protosyn-api/common/introduction.md"
            ],
            hide("Internals" => "protosyn-api/internals.md"),
        ]
    ],
)

deploydocs(
    repo = "github.com/sergio-santos-group/ProtoSyn.jl.git"
)