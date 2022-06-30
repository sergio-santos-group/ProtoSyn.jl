using Printf

@testset verbose = true "Peptides | Mutators $(repeat("-", 37))" begin

    @testset verbose = true "$(@sprintf "%-52s" "Rotamer")" begin
        pose = copy(backup)

        # Set standard deviation of the target rotamer to be 0.0
        r = ProtoSyn.Peptides.load_dunbrack(Float64)
        chi1 = r["MET"][3.14, 3.14][1][chi"1"]
        r["MET"][3.14, 3.14][1][chi"1"] = (chi1[1], 0.0)
        chi2 = r["MET"][3.14, 3.14][1][chi"2"]
        r["MET"][3.14, 3.14][1][chi"2"] = (chi2[1], 0.0)
        chi3 = r["MET"][3.14, 3.14][1][chi"3"]
        r["MET"][3.14, 3.14][1][chi"3"] = (chi3[1], 0.0)
        
        m = ProtoSyn.Peptides.Mutators.RotamerMutator(
            r, 1.0, 1, an"CA" & rid"2", false)

        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CG"]) ≈ 3.141387 atol = 1e-5
        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["SD"]) ≈ -3.141542 atol = 1e-5
        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CE"]) ≈ -3.14153 atol = 1e-5
        @test pose.state.i2c == false
        m(pose)
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

    @testset verbose = true "$(@sprintf "%-52s" "Design")" begin
        pose = copy(backup)
        
        # Force mutation to ALA
        m = ProtoSyn.Peptides.Mutators.DesignMutator(1.0, res_lib, an"CA" & rid"2")
        for aa in keys(ProtoSyn.Peptides.available_aminoacids)
            m.searchable_aminoacids[aa] = false
        end
        m.searchable_aminoacids['A'] = true

        @test pose.graph[1][2].name == "MET"
        @test pose.graph[1][2].size == 17
        @test length(pose.graph[1][2].items) == 17
        @test pose.graph[1][2]["HB3"] === nothing
        @test pose.graph[1][2]["CG"] !== nothing
        @test pose.state.i2c == false
        m(pose)
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test pose.graph[1][2].name == "ALA"
        @test pose.graph[1][2].size == 10
        @test length(pose.graph[1][2].items) == 10
        @test pose.graph[1][2]["HB3"] !== nothing
        @test pose.graph[1][2]["CG"] === nothing
    end
end