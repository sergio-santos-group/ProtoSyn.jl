@testset verbose = true "Peptides | Rotamers $(repeat("-", 37))" begin
    @testset verbose = true "$(@sprintf "%-54s" "Load Dunbrack")" begin
        rot_lib = ProtoSyn.Peptides.load_dunbrack(Float64)
        @test length(keys(rot_lib)) == 19
        @test rot_lib.count == 19
        @test isa(rot_lib["VAL"], ProtoSyn.Peptides.BBD_RotamerLibrary)
        @test isa(rot_lib["VAL"][0.1, 0.1], ProtoSyn.Peptides.BBI_RotamerLibrary)
        @test length(rot_lib["VAL"][0.1, 0.1].rotamers) == 3
        @test length(rot_lib["VAL"][0.1, 0.1].weights) == 3
        rot = rot_lib["VAL"][0.1, 0.1][1]
        @test isa(rot, ProtoSyn.Peptides.Rotamer)
        @test length(values(rot.chis)) == 1
        @test rot.name == "VAL"
    end

    @testset verbose = true "$(@sprintf "%-54s" "Applying Rotamer")" begin
        pose = copy(backup)
        rot = ProtoSyn.Peptides.Rotamer("MET",
            Dict{AbstractSelection, Tuple{Union{Nothing, Float64}, Float64}}(
                chi"1" => (deg2rad(63.9), deg2rad(0.0)),
                chi"2" => (deg2rad(-172.6), deg2rad(0.0)),
                chi"3" => (deg2rad(72.0), deg2rad(0.0))
            )
        )

        @test rad2deg(ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CG"])) ≈ 179.98821691726314 atol = 1e-5
        @test rad2deg(ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["SD"])) ≈ -179.99709776308765 atol = 1e-5
        @test rad2deg(ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CE"])) ≈ -179.9964102137335 atol = 1e-5
        ProtoSyn.Peptides.apply!(pose.state, rot, pose.graph[1][2])
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test rad2deg(ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CG"])) ≈ 63.9 atol = 1e-5
        @test rad2deg(ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["SD"])) ≈ -172.6 atol = 1e-5
        @test rad2deg(ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CE"])) ≈ 72.0 atol = 1e-5
    end

    @testset verbose = true "$(@sprintf "%-54s" "Get Rotamer")" begin
        pose = copy(backup)
        rot = ProtoSyn.Peptides.get_rotamer(pose, pose.graph[1][2])
        @test rot[chi"1"][1] ≈ 3.141387 atol = 1e-5
        @test rot[chi"2"][1] ≈ -3.14154 atol = 1e-5
        @test rot[chi"3"][1] ≈ -3.14153 atol = 1e-5
        @test rot[chi"1"][2] == 0.0
        @test rot[chi"2"][2] == 0.0
        @test rot[chi"3"][2] == 0.0
    end
end