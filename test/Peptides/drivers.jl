using Printf

@testset verbose = true "Peptides | Drivers $(repeat("-", 35))" begin

    @testset verbose = true "$(@sprintf "%-52s" "Rotamer Blitz")" begin
        pose = copy(backup)
        energy_function = ProtoSyn.Calculators.EnergyFunction([
            ProtoSyn.Calculators.TorchANI.get_default_torchani_ensemble()
        ])

        # Set standard deviation of the target rotamer to be 0.0
        r = copy(rot_lib)
        chi1 = r["MET"][3.14, 3.14][1].chis[ProtoSyn.Peptides.Dihedral.Chi1()]
        r["MET"][3.14, 3.14][1].chis[ProtoSyn.Peptides.Dihedral.Chi1()] = (chi1[1], 0.0)
        chi2 = r["MET"][3.14, 3.14][1].chis[ProtoSyn.Peptides.Dihedral.Chi2()]
        r["MET"][3.14, 3.14][1].chis[ProtoSyn.Peptides.Dihedral.Chi2()] = (chi2[1], 0.0)
        chi3 = r["MET"][3.14, 3.14][1].chis[ProtoSyn.Peptides.Dihedral.Chi3()]
        r["MET"][3.14, 3.14][1].chis[ProtoSyn.Peptides.Dihedral.Chi3()] = (chi3[1], 0.0)
        
        rb = ProtoSyn.Peptides.Drivers.RotamerBlitz(energy_function, r, 1, 1, nothing, rid"2")

        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CG"]) ≈ 3.141387 atol = 1e-5
        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["SD"]) ≈ -3.141542 atol = 1e-5
        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CE"]) ≈ -3.14153 atol = 1e-5
        @test pose.state.i2c == false
        rb(pose)
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CG"]) ≈ 1.1152653920243765 atol = 1e-5
        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["SD"]) ≈ -3.0124382889422128 atol = 1e-5
        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CE"]) ≈ 1.2566370614359168 atol = 1e-5
        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CG"]) ≈ chi1[1] atol = 1e-5
        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["SD"]) ≈ chi2[1] atol = 1e-5
        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CE"]) ≈ chi3[1] atol = 1e-5
    end
end