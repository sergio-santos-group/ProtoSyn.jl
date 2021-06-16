using Printf

@testset verbose = true "Peptides | Calculators $(repeat("-", 31))" begin

    @testset verbose = true "$(@sprintf "%-52s" "Sidechain Clash Restraint")" begin
        pose = copy(backup)

        # With dynamic mask
        c1 = ProtoSyn.Peptides.Calculators.Restraints.get_default_sidechain_clash_restraint()
        r = c1.calc(pose, false; c1.settings...)
        @test r[1] == 0.0
        @test sum(r[2]) == 0.0

        # With fixed mask
        mask = ProtoSyn.Calculators.intra_residue_mask(pose, !an"^C$|^N$|^CA$|^H$|^O$"r)
        c2 = ProtoSyn.Peptides.Calculators.Restraints.get_default_sidechain_clash_restraint(mask = mask)
        c2.settings[:d1] = 10.0
        c2.settings[:d2] = 15.0
        r = c2.calc(pose, true; c2.settings...)
        @test r[1] ≈ 17761.286101800346 atol = 1e-5
        @test sum(r[2]) ≈ 2.2737367544323206e-13 atol = 1e-10
    end

    @testset verbose = true "$(@sprintf "%-52s" "Contact Restraint")" begin
        pose = copy(backup)

        # With fixed map
        c1 = ProtoSyn.Peptides.Calculators.Restraints.get_default_contact_restraint("../Core/cmap-test-1.txt")
        r = c1.calc(pose, false; c1.settings...)
        @test r[1] == 0.0
        @test sum(r[2]) == 0.0
    end

    @testset verbose = true "$(@sprintf "%-52s" "Ca-Ca Clash Restraint")" begin
        pose = copy(backup)

        # With dynamic mask
        c1 = ProtoSyn.Peptides.Calculators.Restraints.get_default_ca_clash_restraint()
        r = c1.calc(pose, false; c1.settings...)
        @test r[1] == 0.0
        @test sum(r[2]) == 0.0

        # With fixed mask
        mask = ProtoSyn.Calculators.diagonal_mask(pose, an"CA")
        c2 = ProtoSyn.Peptides.Calculators.Restraints.get_default_ca_clash_restraint(mask = mask)
        c2.settings[:d1] = 10.0
        c2.settings[:d2] = 15.0
        r = c2.calc(pose, true; c2.settings...)
        @test r[1] ≈ 459.91112063176973 atol = 1e-5
        @test sum(r[2]) ≈ -7.105427357601002e-15 atol = 1e-10
    end

    @testset verbose = true "$(@sprintf "%-52s" "Caterpillar Solvation")" begin
        pose = copy(backup)

        # With dynamic mask
        c1 = ProtoSyn.Peptides.Calculators.Caterpillar.get_default_solvation_energy()
        r = c1.calc(pose, false; c1.settings...)
        @test r[1] ≈ 104.99999999992491 atol = 1e-5
        @test r[2] === nothing
    end
end