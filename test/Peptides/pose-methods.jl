@testset verbose = true "Pose methods" begin
    @testset verbose = true "Append fragment (from fragment, to end)" begin
        pose = copy(backup)
        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.state.size == 39
        @test frag.state.size == 30
        @test length(collect(eachresidue(pose.graph))) == 3
        @test length(collect(eachresidue(frag.graph))) == 3

        ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][end], res_lib, frag)
        @test pose.state.size == (39 + 30)
        @test length(collect(eachresidue(pose.graph))) == (3 + 3)
        @test pose.graph[1][4].parent == pose.graph[1][3]
        @test pose.graph[1][4]["N"].parent == pose.graph[1][3]["C"]
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test collect(pose.state[pose.graph[1][4]["N"]].t) ≈ [9.540142155877003, -6.733599261965266, -2.6094253145331535e-6]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.graph[1][4]["N"].id == 40
        @test pose.graph[1][4].id == 4
        @test pose.graph[1][4]["N"].ascendents == (40, 38, 27, 25)
    end

    @testset verbose = true "Append fragment (from fragment, to unbond)" begin
        pose = copy(backup)

        ProtoSyn.unbond!(pose, pose.graph[1][2]["C"], pose.graph[1][3]["N"])

        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.state.size == 39
        @test frag.state.size == 30
        @test length(collect(eachresidue(pose.graph))) == 3
        @test length(collect(eachresidue(frag.graph))) == 3

        ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][2], res_lib, frag)
        @test pose.graph[1][3].name == "ALA"
        @test pose.graph[1][4].name == "ALA"
        @test pose.graph[1][5].name == "ALA"
        @test pose.state.size == (39 + 30)
        @test length(collect(eachresidue(pose.graph))) == (3 + 3)
        @test pose.graph[1][3].parent == pose.graph[1][2]
        @test pose.graph[1][3]["N"].parent == pose.graph[1][2]["C"]
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test collect(pose.state[pose.graph[1][3]["N"]].t) ≈ [6.377553820615265, -5.088979503045568, -2.1339193688630605e-6]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898658, -1.0392305459771767, -3.6000006795018745e-7]
        @test pose.graph[1][3]["N"].id == 25
        @test pose.graph[1][3].id == 3
        @test pose.graph[1][3]["N"].ascendents == (25, 23, 10, 8)

        @test pose.graph[1][6]["N"].parent == ProtoSyn.root(pose.graph)
        @test pose.graph[1][6].parent == ProtoSyn.root(pose.graph).container
    end

    @testset verbose = true "Insert fragment (from fragment, to middle)" begin
        pose = copy(backup)
        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.state.size == 39
        @test frag.state.size == 30
        @test length(collect(eachresidue(pose.graph))) == 3
        @test length(collect(eachresidue(frag.graph))) == 3
        @test pose.graph[1][2].name == "MET"
        @test pose.graph[1][2].parent == pose.graph[1][1]
        @test pose.graph[1][2]["N"].parent == pose.graph[1][1]["C"]

        ProtoSyn.Peptides.insert_fragment!(pose, pose.graph[1][2], res_lib, frag)
        @test pose.graph[1][2].name == "ALA"
        @test pose.graph[1][3].name == "ALA"
        @test pose.graph[1][4].name == "ALA"
        @test pose.state.size == (39 + 30)
        @test length(collect(eachresidue(pose.graph))) == (3 + 3)
        @test pose.graph[1][2].parent == pose.graph[1][1]
        @test pose.graph[1][2]["N"].parent == pose.graph[1][1]["C"]
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test collect(pose.state[pose.graph[1][2]["N"]].t) ≈ [3.7626151211853136, -2.6664117066805653, -8.203956751375917e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898658, -1.0392305459771767, -3.6000006795018745e-7]
        @test pose.graph[1][2]["N"].id == 8
        @test pose.graph[1][2].id == 2
        @test pose.graph[1][4]["N"].ascendents == (28, 26, 20, 18)

        @test pose.graph[1][5].name == "MET"
        @test pose.graph[1][5].parent == pose.graph[1][4]
        @test pose.graph[1][5]["N"].parent == pose.graph[1][4]["C"]
        @test collect(pose.state[pose.graph[1][5]["N"]].t) ≈ [12.155080855306954, -9.156167058330269, -3.922949007961368e-6]
    end

    @testset verbose = true "Insert fragment (from fragment, to start)" begin
        pose = copy(backup)
        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.state.size == 39
        @test frag.state.size == 30
        @test length(collect(eachresidue(pose.graph))) == 3
        @test length(collect(eachresidue(frag.graph))) == 3
        @test pose.graph[1][1].name == "GLY"
        @test pose.graph[1][1].parent == ProtoSyn.root(pose.graph).container
        @test pose.graph[1][1]["N"].parent == ProtoSyn.root(pose.graph)

        ProtoSyn.Peptides.insert_fragment!(pose, pose.graph[1][1], res_lib, frag)
        @test pose.graph[1][1].name == "ALA"
        @test pose.graph[1][2].name == "ALA"
        @test pose.graph[1][3].name == "ALA"
        @test pose.state.size == (39 + 30)
        @test length(collect(eachresidue(pose.graph))) == (3 + 3)
        @test pose.graph[1][1].parent == ProtoSyn.root(pose.graph).container
        @test pose.graph[1][1]["N"].parent == ProtoSyn.root(pose.graph)
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test collect(pose.state[pose.graph[1][2]["N"]].t) ≈ [3.781926783267048, -2.646115855671686, -7.970254511993344e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.graph[1][2]["N"].id == 11
        @test pose.graph[1][2].id == 2
        @test pose.graph[1][4]["N"].ascendents == (31, 29, 23, 21)

        @test pose.graph[1][4].name == "GLY"
        @test pose.graph[1][4].parent == pose.graph[1][3]
        @test pose.graph[1][4]["N"].parent == pose.graph[1][3]["C"]
        @test collect(pose.state[pose.graph[1][5]["N"]].t) ≈ [12.262986826440596, -9.010268073403653, -3.764224971305979e-6]
    end

    # --------------------------------------------------------------------------

    @testset verbose = true "Append fragment (from derivation, to end)" begin
        pose = copy(backup)
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.state.size == 39
        @test frag.state.size == 30
        @test length(collect(eachresidue(pose.graph))) == 3
        @test length(collect(eachresidue(frag.graph))) == 3

        ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][end], res_lib, seq"AAA")
        @test pose.state.size == (39 + 30)
        @test length(collect(eachresidue(pose.graph))) == (3 + 3)
        @test pose.graph[1][4].parent == pose.graph[1][3]
        @test pose.graph[1][4]["N"].parent == pose.graph[1][3]["C"]
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test collect(pose.state[pose.graph[1][4]["N"]].t) ≈ [9.540142155877003, -6.733599261965266, -2.6094253145331535e-6]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.graph[1][4]["N"].id == 40
        @test pose.graph[1][4].id == 4
        @test pose.graph[1][4]["N"].ascendents == (40, 38, 27, 25)
    end

    @testset verbose = true "Append fragment (from derivation, to unbond)" begin
        pose = copy(backup)

        ProtoSyn.unbond!(pose, pose.graph[1][2]["C"], pose.graph[1][3]["N"])

        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.state.size == 39
        @test frag.state.size == 30
        @test length(collect(eachresidue(pose.graph))) == 3
        @test length(collect(eachresidue(frag.graph))) == 3

        ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][2], res_lib, seq"AAA")
        @test pose.graph[1][3].name == "ALA"
        @test pose.graph[1][4].name == "ALA"
        @test pose.graph[1][5].name == "ALA"
        @test pose.state.size == (39 + 30)
        @test length(collect(eachresidue(pose.graph))) == (3 + 3)
        @test pose.graph[1][3].parent == pose.graph[1][2]
        @test pose.graph[1][3]["N"].parent == pose.graph[1][2]["C"]
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test collect(pose.state[pose.graph[1][3]["N"]].t) ≈ [6.377553820615265, -5.088979503045568, -2.1339193688630605e-6]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898658, -1.0392305459771767, -3.6000006795018745e-7]
        @test pose.graph[1][3]["N"].id == 25
        @test pose.graph[1][3].id == 3
        @test pose.graph[1][3]["N"].ascendents == (25, 23, 10, 8)

        @test pose.graph[1][6]["N"].parent == ProtoSyn.root(pose.graph)
        @test pose.graph[1][6].parent == ProtoSyn.root(pose.graph).container
    end

    @testset verbose = true "Insert fragment (from derivation, to middle)" begin
        pose = copy(backup)
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.state.size == 39
        @test frag.state.size == 30
        @test length(collect(eachresidue(pose.graph))) == 3
        @test length(collect(eachresidue(frag.graph))) == 3
        @test pose.graph[1][2].name == "MET"
        @test pose.graph[1][2].parent == pose.graph[1][1]
        @test pose.graph[1][2]["N"].parent == pose.graph[1][1]["C"]

        ProtoSyn.Peptides.insert_fragment!(pose, pose.graph[1][2], res_lib, seq"AAA")
        @test pose.graph[1][2].name == "ALA"
        @test pose.graph[1][3].name == "ALA"
        @test pose.graph[1][4].name == "ALA"
        @test pose.state.size == (39 + 30)
        @test length(collect(eachresidue(pose.graph))) == (3 + 3)
        @test pose.graph[1][2].parent == pose.graph[1][1]
        @test pose.graph[1][2]["N"].parent == pose.graph[1][1]["C"]
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test collect(pose.state[pose.graph[1][2]["N"]].t) ≈ [3.7626151211853136, -2.6664117066805653, -8.203956751375917e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898658, -1.0392305459771767, -3.6000006795018745e-7]
        @test pose.graph[1][2]["N"].id == 8
        @test pose.graph[1][2].id == 2
        @test pose.graph[1][4]["N"].ascendents == (28, 26, 20, 18)

        @test pose.graph[1][5].name == "MET"
        @test pose.graph[1][5].parent == pose.graph[1][4]
        @test pose.graph[1][5]["N"].parent == pose.graph[1][4]["C"]
        @test collect(pose.state[pose.graph[1][5]["N"]].t) ≈ [12.155080855306954, -9.156167058330269, -3.922949007961368e-6]
    end

    @testset verbose = true "Insert fragment (from derivation, to start)" begin
        pose = copy(backup)
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.state.size == 39
        @test frag.state.size == 30
        @test length(collect(eachresidue(pose.graph))) == 3
        @test length(collect(eachresidue(frag.graph))) == 3
        @test pose.graph[1][1].name == "GLY"
        @test pose.graph[1][1].parent == ProtoSyn.root(pose.graph).container
        @test pose.graph[1][1]["N"].parent == ProtoSyn.root(pose.graph)

        ProtoSyn.Peptides.insert_fragment!(pose, pose.graph[1][1], res_lib, seq"AAA")
        @test pose.graph[1][1].name == "ALA"
        @test pose.graph[1][2].name == "ALA"
        @test pose.graph[1][3].name == "ALA"
        @test pose.state.size == (39 + 30)
        @test length(collect(eachresidue(pose.graph))) == (3 + 3)
        @test pose.graph[1][1].parent == ProtoSyn.root(pose.graph).container
        @test pose.graph[1][1]["N"].parent == ProtoSyn.root(pose.graph)
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test collect(pose.state[pose.graph[1][2]["N"]].t) ≈ [3.781926783267048, -2.646115855671686, -7.970254511993344e-7]
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7]
        @test pose.graph[1][2]["N"].id == 11
        @test pose.graph[1][2].id == 2
        @test pose.graph[1][4]["N"].ascendents == (31, 29, 23, 21)

        @test pose.graph[1][4].name == "GLY"
        @test pose.graph[1][4].parent == pose.graph[1][3]
        @test pose.graph[1][4]["N"].parent == pose.graph[1][3]["C"]
        @test collect(pose.state[pose.graph[1][5]["N"]].t) ≈ [12.262986826440596, -9.010268073403653, -3.764224971305979e-6]
    end
end