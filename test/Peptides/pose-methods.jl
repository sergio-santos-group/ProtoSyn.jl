using Printf

@testset verbose = true "Peptides | Pose methods $(repeat("-", 33))" begin

    @testset verbose = true "$(@sprintf "%-54s" "Identify terminals")" begin
        pose = copy(backup)

        @test ProtoSyn.Peptides.is_C_terminal(pose.graph[1, 1]) === false
        @test ProtoSyn.Peptides.is_C_terminal(pose.graph[1, 3]) === true
        @test ProtoSyn.Peptides.is_N_terminal(pose.graph[1, 1]) === true
        @test ProtoSyn.Peptides.is_N_terminal(pose.graph[1, 3]) === false
        @test ProtoSyn.Peptides.identify_c_terminal(pose.graph[1]) === pose.graph[1, 3, "C"]
    end

    @testset verbose = true "$(@sprintf "%-54s" "Append fragment (from fragment, to end)")" begin
        pose = copy(backup)
        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
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
        @test collect(pose.state[pose.graph[1][4]["N"]].t) ≈ [9.540142155877003, -6.733599261965266, -2.6094253145331535e-6] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
        @test pose.graph[1][4]["N"].id == 40
        @test pose.graph[1][4].id == 4
        @test pose.graph[1][4]["N"].ascendents == (40, 38, 27, 25)
    end

    @testset verbose = true "$(@sprintf "%-54s" "Append fragment (from fragment, to unbond)")" begin
        pose = copy(backup)

        ProtoSyn.unbond!(pose, pose.graph[1][2]["C"], pose.graph[1][3]["N"])

        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
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
        @test collect(pose.state[pose.graph[1][3]["N"]].t) ≈ [6.377553820615265, -5.088979503045568, -2.1339193688630605e-6] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898658, -1.0392305459771767, -3.6000006795018745e-7] atol = 1e-5
        @test pose.graph[1][3]["N"].id == 25
        @test pose.graph[1][3].id == 3
        @test pose.graph[1][3]["N"].ascendents == (25, 23, 10, 8)

        @test pose.graph[1][6]["N"].parent == ProtoSyn.root(pose.graph)
        @test pose.graph[1][6].parent == ProtoSyn.root(pose.graph).container
    end

    @testset verbose = true "$(@sprintf "%-54s" "Insert fragment (from fragment, to middle)")" begin
        pose = copy(backup)
        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
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
        @test collect(pose.state[pose.graph[1][2]["N"]].t) ≈ [3.7626151211853136, -2.6664117066805653, -8.203956751375917e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898658, -1.0392305459771767, -3.6000006795018745e-7] atol = 1e-5
        @test pose.graph[1][2]["N"].id == 8
        @test pose.graph[1][2].id == 2
        @test pose.graph[1][4]["N"].ascendents == (28, 26, 20, 18)

        @test pose.graph[1][5].name == "MET"
        @test pose.graph[1][5].parent == pose.graph[1][4]
        @test pose.graph[1][5]["N"].parent == pose.graph[1][4]["C"]
        @test collect(pose.state[pose.graph[1][5]["N"]].t) ≈ [12.155080855306954, -9.156167058330269, -3.922949007961368e-6] atol = 1e-5
    end

    @testset verbose = true "$(@sprintf "%-54s" "Insert fragment (from fragment, to start)")" begin
        pose = copy(backup)
        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
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
        @test collect(pose.state[pose.graph[1][2]["N"]].t) ≈ [3.781926783267048, -2.646115855671686, -7.970254511993344e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
        @test pose.graph[1][2]["N"].id == 11
        @test pose.graph[1][2].id == 2
        @test pose.graph[1][4]["N"].ascendents == (31, 29, 23, 21)

        @test pose.graph[1][4].name == "GLY"
        @test pose.graph[1][4].parent == pose.graph[1][3]
        @test pose.graph[1][4]["N"].parent == pose.graph[1][3]["C"]
        @test collect(pose.state[pose.graph[1][5]["N"]].t) ≈ [12.262986826440596, -9.010268073403653, -3.764224971305979e-6] atol = 1e-5
    end

    # --------------------------------------------------------------------------

    @testset verbose = true "$(@sprintf "%-54s" "Append fragment (from derivation, to end)")" begin
        pose = copy(backup)
        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
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
        @test collect(pose.state[pose.graph[1][4]["N"]].t) ≈ [9.540142155877003, -6.733599261965266, -2.6094253145331535e-6] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
        @test pose.graph[1][4]["N"].id == 40
        @test pose.graph[1][4].id == 4
        @test pose.graph[1][4]["N"].ascendents == (40, 38, 27, 25)
    end

    @testset verbose = true "$(@sprintf "%-54s" "Append fragment (from derivation, to unbond)")" begin
        pose = copy(backup)
        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        ProtoSyn.unbond!(pose, pose.graph[1][2]["C"], pose.graph[1][3]["N"])

        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
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
        @test collect(pose.state[pose.graph[1][3]["N"]].t) ≈ [6.377553820615265, -5.088979503045568, -2.1339193688630605e-6] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898658, -1.0392305459771767, -3.6000006795018745e-7] atol = 1e-5
        @test pose.graph[1][3]["N"].id == 25
        @test pose.graph[1][3].id == 3
        @test pose.graph[1][3]["N"].ascendents == (25, 23, 10, 8)

        @test pose.graph[1][6]["N"].parent == ProtoSyn.root(pose.graph)
        @test pose.graph[1][6].parent == ProtoSyn.root(pose.graph).container
    end

    @testset verbose = true "$(@sprintf "%-54s" "Insert fragment (from derivation, to middle)")" begin
        pose = copy(backup)
        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
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
        @test collect(pose.state[pose.graph[1][2]["N"]].t) ≈ [3.7626151211853136, -2.6664117066805653, -8.203956751375917e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898658, -1.0392305459771767, -3.6000006795018745e-7] atol = 1e-5
        @test pose.graph[1][2]["N"].id == 8
        @test pose.graph[1][2].id == 2
        @test pose.graph[1][4]["N"].ascendents == (28, 26, 20, 18)

        @test pose.graph[1][5].name == "MET"
        @test pose.graph[1][5].parent == pose.graph[1][4]
        @test pose.graph[1][5]["N"].parent == pose.graph[1][4]["C"]
        @test collect(pose.state[pose.graph[1][5]["N"]].t) ≈ [12.155080855306954, -9.156167058330269, -3.922949007961368e-6] atol = 1e-5
    end

    @testset verbose = true "$(@sprintf "%-54s" "Insert fragment (from derivation, to start)")" begin
        pose = copy(backup)
        frag = ProtoSyn.fragment(res_lib, seq"AAA")
        @test collect(frag.state[frag.graph[1]["N"]].t) ≈ [0.5999998935898667, -1.0392305459771758, -3.6000006815716016e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
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
        @test collect(pose.state[pose.graph[1][2]["N"]].t) ≈ [3.781926783267048, -2.646115855671686, -7.970254511993344e-7] atol = 1e-5
        @test collect(pose.state[pose.graph[1][1]["N"]].t) ≈ [0.5999998935898657, -1.0392305459771765, -3.600000679501874e-7] atol = 1e-5
        @test pose.graph[1][2]["N"].id == 11
        @test pose.graph[1][2].id == 2
        @test pose.graph[1][4]["N"].ascendents == (31, 29, 23, 21)

        @test pose.graph[1][4].name == "GLY"
        @test pose.graph[1][4].parent == pose.graph[1][3]
        @test pose.graph[1][4]["N"].parent == pose.graph[1][3]["C"]
        @test collect(pose.state[pose.graph[1][5]["N"]].t) ≈ [12.262986826440596, -9.010268073403653, -3.764224971305979e-6] atol = 1e-5
    end

    # --------------------------------------------------------------------------
    # Mutation

    @testset verbose = true "$(@sprintf "%-54s" "Mutate")" begin
        pose = copy(backup)
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 39
        @test length(pose.graph[1].items) == 3
        @test collect(pose.state[pose.graph[1][2]["CB"]].t) ≈ [3.7922803002732848, -4.787107925441947, -1.2320027833649587] atol = 1e-5
        @test pose.graph[1][3] in pose.graph[1][2].children
        @test pose.graph[1][3].parent == pose.graph[1][2]
        @test pose.graph[1][3]["N"] in pose.graph[1][2]["C"].children
        @test pose.graph[1][3]["N"].parent == pose.graph[1][2]["C"]
        @test pose.graph[1][2] in pose.graph[1][1].children
        @test pose.graph[1][2].parent == pose.graph[1][1]
        @test pose.graph[1][2]["N"] in pose.graph[1][1]["C"].children
        @test pose.graph[1][2]["N"].parent == pose.graph[1][1]["C"]
        @test pose.graph[1][3]["N"].ascendents == (25, 23, 10, 8)

        ProtoSyn.Peptides.mutate!(pose, pose.graph[1][2], res_lib, seq"I")

        @test pose.state.i2c == true
        @test pose.graph[1][2].name == "ILE"
        @test pose.state.size == 41
        @test length(pose.graph[1].items) == 3
        @test pose.graph[1][3] in pose.graph[1][2].children
        @test pose.graph[1][3].parent == pose.graph[1][2]
        @test pose.graph[1][3]["N"] in pose.graph[1][2]["C"].children
        @test pose.graph[1][3]["N"].parent == pose.graph[1][2]["C"]
        @test pose.graph[1][2] in pose.graph[1][1].children
        @test pose.graph[1][2].parent == pose.graph[1][1]
        @test pose.graph[1][2]["N"] in pose.graph[1][1]["C"].children
        @test pose.graph[1][2]["N"].parent == pose.graph[1][1]["C"]
        @test pose.graph[1][3]["N"].ascendents == (27, 25, 10, 8)

        sync!(pose)
        @test collect(pose.state[pose.graph[1][2]["CB"]].t) ≈ [3.770863352158185, -4.75167702051048, -1.2450030977749311] atol = 1e-5
    end

    @testset verbose = true "$(@sprintf "%-54s" "Force mutate")" begin
        pose = copy(backup)
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 39
        @test length(pose.graph[1].items) == 3
        @test collect(pose.state[pose.graph[1][2]["CB"]].t) ≈ [3.7922803002732848, -4.787107925441947, -1.2320027833649587] atol = 1e-5
        @test pose.graph[1][3] in pose.graph[1][2].children
        @test pose.graph[1][3].parent == pose.graph[1][2]
        @test pose.graph[1][3]["N"] in pose.graph[1][2]["C"].children
        @test pose.graph[1][3]["N"].parent == pose.graph[1][2]["C"]
        @test pose.graph[1][2] in pose.graph[1][1].children
        @test pose.graph[1][2].parent == pose.graph[1][1]
        @test pose.graph[1][2]["N"] in pose.graph[1][1]["C"].children
        @test pose.graph[1][2]["N"].parent == pose.graph[1][1]["C"]
        @test pose.graph[1][3]["N"].ascendents == (25, 23, 10, 8)

        ProtoSyn.Peptides.force_mutate!(pose, pose.graph[1][2], res_lib, seq"I")

        @test pose.state.i2c == true
        @test pose.graph[1][2].name == "ILE"
        @test pose.state.size == 41
        @test length(pose.graph[1].items) == 3
        @test pose.graph[1][3] in pose.graph[1][2].children
        @test pose.graph[1][3].parent == pose.graph[1][2]
        @test pose.graph[1][3]["N"] in pose.graph[1][2]["C"].children
        @test pose.graph[1][3]["N"].parent == pose.graph[1][2]["C"]
        @test pose.graph[1][2] in pose.graph[1][1].children
        @test pose.graph[1][2].parent == pose.graph[1][1]
        @test pose.graph[1][2]["N"] in pose.graph[1][1]["C"].children
        @test pose.graph[1][2]["N"].parent == pose.graph[1][1]["C"]
        @test pose.graph[1][3]["N"].ascendents == (27, 25, 10, 8)

        sync!(pose)
        @test collect(pose.state[pose.graph[1][2]["CB"]].t) ≈ [3.770863352158185, -4.75167702051048, -1.2450030977749311] atol = 1e-5
    end

    # --------------------------------------------------------------------------
    # Sidechain removal /add

    @testset verbose = true "$(@sprintf "%-54s" "Remove sidechain (no selection)")" begin
        pose = copy(backup)
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 39
        @test pose.graph[1][2]["CB"] in pose.graph[1][2]["CA"].children
        @test length(pose.graph[1][2]["CA"].children) === 3
        @test pose.graph[1][2]["HA1"] === nothing
        @test pose.graph[1][2]["HA2"] === nothing

        ProtoSyn.Peptides.remove_sidechains!(pose, res_lib)
        @test pose.state.i2c == false
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 21
        @test pose.graph[1][2]["CB"] === nothing
        @test length(pose.graph[1][2]["CA"].children) === 3
        @test pose.graph[1][2]["HA1"] in pose.graph[1][2]["CA"].children
        @test pose.graph[1][2]["HA2"] in pose.graph[1][2]["CA"].children
        @test collect(pose.state[pose.graph[1][2]["HA1"]].t) ≈ [3.9132885320002266, -4.544421783453538, 0.8899977141815026] atol = 1e-5
    end

    @testset verbose = true "$(@sprintf "%-54s" "Remove sidechain (with selection)")" begin
        pose = copy(backup)
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 39
        @test pose.graph[1][2]["CB"] in pose.graph[1][2]["CA"].children
        @test length(pose.graph[1][2]["CA"].children) === 3
        @test pose.graph[1][2]["HA1"] === nothing
        @test pose.graph[1][2]["HA2"] === nothing
        @test pose.graph[1][3].name == "GLU"
        @test pose.graph[1][3]["CB"] in pose.graph[1][3]["CA"].children
        @test length(pose.graph[1][3]["CA"].children) === 3
        @test pose.graph[1][3]["HA1"] === nothing
        @test pose.graph[1][3]["HA2"] === nothing

        ProtoSyn.Peptides.remove_sidechains!(pose, res_lib, rid"2")
        @test pose.state.i2c == false
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 29
        @test pose.graph[1][2]["CB"] === nothing
        @test length(pose.graph[1][2]["CA"].children) === 3
        @test pose.graph[1][2]["HA1"] in pose.graph[1][2]["CA"].children
        @test pose.graph[1][2]["HA2"] in pose.graph[1][2]["CA"].children
        @test collect(pose.state[pose.graph[1][2]["HA1"]].t) ≈ [3.9132885320002266, -4.544421783453538, 0.8899977141815026] atol = 1e-5
        @test pose.graph[1][3].name == "GLU"
        @test pose.graph[1][3]["CB"] in pose.graph[1][3]["CA"].children
        @test length(pose.graph[1][3]["CA"].children) === 3
        @test pose.graph[1][3]["HA1"] === nothing
        @test pose.graph[1][3]["HA2"] === nothing
    end

    @testset verbose = true "$(@sprintf "%-54s" "Force remove sidechain (no selection)")" begin
        pose = copy(backup)
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 39
        @test pose.graph[1][2]["CB"] in pose.graph[1][2]["CA"].children
        @test length(pose.graph[1][2]["CA"].children) === 3
        @test pose.graph[1][2]["HA1"] === nothing
        @test pose.graph[1][2]["HA2"] === nothing

        ProtoSyn.Peptides.force_remove_sidechains!(pose)
        @test pose.state.i2c == false
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 15
        @test pose.graph[1][2]["CB"] === nothing
        @test length(pose.graph[1][2]["CA"].children) === 1
        @test pose.graph[1][2]["HA1"] === nothing
        @test pose.graph[1][2]["HA2"] === nothing
    end

    @testset verbose = true "$(@sprintf "%-54s" "Force remove sidechain (with selection)")" begin
        pose = copy(backup)
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 39
        @test pose.graph[1][2]["CB"] in pose.graph[1][2]["CA"].children
        @test length(pose.graph[1][2]["CA"].children) === 3
        @test pose.graph[1][2]["HA1"] === nothing
        @test pose.graph[1][2]["HA2"] === nothing
        @test pose.graph[1][3].name == "GLU"
        @test pose.graph[1][3]["CB"] in pose.graph[1][3]["CA"].children
        @test length(pose.graph[1][3]["CA"].children) === 3
        @test pose.graph[1][3]["HA1"] === nothing
        @test pose.graph[1][3]["HA2"] === nothing

        ProtoSyn.Peptides.force_remove_sidechains!(pose, rid"2")
        @test pose.state.i2c == false
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 27
        @test pose.graph[1][2]["CB"] === nothing
        @test length(pose.graph[1][2]["CA"].children) === 1
        @test pose.graph[1][2]["HA1"] === nothing
        @test pose.graph[1][2]["HA2"] === nothing
        @test pose.graph[1][3].name == "GLU"
        @test pose.graph[1][3]["CB"] in pose.graph[1][3]["CA"].children
        @test length(pose.graph[1][3]["CA"].children) === 3
        @test pose.graph[1][3]["HA1"] === nothing
        @test pose.graph[1][3]["HA2"] === nothing
    end

    @testset verbose = true "$(@sprintf "%-54s" "Add sidechain (no selection)")" begin
        pose = copy(backup)
        ProtoSyn.Peptides.remove_sidechains!(pose, res_lib)
        @test pose.state.i2c == false
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 21
        @test pose.graph[1][2]["CB"] === nothing
        @test length(pose.graph[1][2]["CA"].children) === 3
        @test pose.graph[1][2]["HA1"] in pose.graph[1][2]["CA"].children
        @test pose.graph[1][2]["HA2"] in pose.graph[1][2]["CA"].children
        @test collect(pose.state[pose.graph[1][2]["HA1"]].t) ≈ [3.9132885320002266, -4.544421783453538, 0.8899977141815026] atol = 1e-5
    
        ProtoSyn.Peptides.add_sidechains!(pose, res_lib)
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 39
        @test pose.graph[1][2]["CB"] in pose.graph[1][2]["CA"].children
        @test length(pose.graph[1][2]["CA"].children) === 3
        @test pose.graph[1][2]["HA1"] === nothing
        @test pose.graph[1][2]["HA2"] === nothing
        @test collect(pose.state[pose.graph[1][2]["CB"]].t) ≈ [3.792279900244243, -4.787108074155822, -1.2320025370097143] atol = 1e-5
    end

    @testset verbose = true "$(@sprintf "%-54s" "Add sidechain (with selection)")" begin
        pose = copy(backup)
        ProtoSyn.Peptides.remove_sidechains!(pose, res_lib)
        @test pose.state.i2c == false
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 21
        @test pose.graph[1][2]["CB"] === nothing
        @test length(pose.graph[1][2]["CA"].children) === 3
        @test pose.graph[1][2]["HA1"] in pose.graph[1][2]["CA"].children
        @test pose.graph[1][2]["HA2"] in pose.graph[1][2]["CA"].children
        @test collect(pose.state[pose.graph[1][2]["HA1"]].t) ≈ [3.9132885320002266, -4.544421783453538, 0.8899977141815026] atol = 1e-5
        @test pose.graph[1][3].name == "GLU"
        @test pose.graph[1][3]["CB"] === nothing
        @test length(pose.graph[1][3]["CA"].children) === 3
        @test pose.graph[1][3]["HA1"] in pose.graph[1][3]["CA"].children
        @test pose.graph[1][3]["HA2"] in pose.graph[1][3]["CA"].children

        ProtoSyn.Peptides.add_sidechains!(pose, res_lib, rid"2")
        @test pose.graph[1][2].name == "MET"
        @test pose.state.size == 31
        @test pose.graph[1][2]["CB"] in pose.graph[1][2]["CA"].children
        @test length(pose.graph[1][2]["CA"].children) === 3
        @test pose.graph[1][2]["HA1"] === nothing
        @test pose.graph[1][2]["HA2"] === nothing
        @test collect(pose.state[pose.graph[1][2]["CB"]].t) ≈ [3.5769993802877447, 3.654001694297366, 1.2320001352913499] atol = 1e-5
        @test pose.graph[1][3].name == "GLU"
        @test pose.graph[1][3]["CB"] === nothing
        @test length(pose.graph[1][3]["CA"].children) === 3
        @test pose.graph[1][3]["HA1"] in pose.graph[1][3]["CA"].children
        @test pose.graph[1][3]["HA2"] in pose.graph[1][3]["CA"].children
    end

    # --------------------------------------------------------------------------
    # Cap/uncap

    @testset verbose = true "$(@sprintf "%-54s" "Uncap (no selection)")" begin
        pose = copy(backup)
        @test pose.state.size == 39
        @test length(pose.graph[1][1]["N"].children) == 2
        @test length(pose.graph[1][end]["C"].children) == 1
        @test pose.graph[1][1]["H"] in pose.graph[1][1]["N"].children
        @test pose.graph[1][1]["H1"] === nothing
        @test pose.graph[1][1]["H2"] === nothing
        @test pose.graph[1][1]["H3"] === nothing
        @test pose.graph[1][end]["O"] in pose.graph[1][end]["C"].children
        @test pose.graph[1][end]["OXT"] === nothing

        ProtoSyn.Peptides.uncap!(pose)
        @test pose.state.size == 37
        @test length(pose.graph[1][1]["N"].children) == 1
        @test length(pose.graph[1][end]["C"].children) == 0
        @test pose.graph[1][1]["H"] === nothing
        @test pose.graph[1][1]["H1"] === nothing
        @test pose.graph[1][1]["H2"] === nothing
        @test pose.graph[1][1]["H3"] === nothing
        @test pose.graph[1][end]["O"] === nothing
        @test pose.graph[1][end]["OXT"] === nothing
    end

    @testset verbose = true "$(@sprintf "%-54s" "Uncap (with selection)")" begin
        pose = copy(backup)
        @test pose.state.size == 39
        @test length(pose.graph[1][1]["N"].children) == 2
        @test length(pose.graph[1][end]["C"].children) == 1
        @test pose.graph[1][1]["H"] in pose.graph[1][1]["N"].children
        @test pose.graph[1][1]["H1"] === nothing
        @test pose.graph[1][1]["H2"] === nothing
        @test pose.graph[1][1]["H3"] === nothing
        @test pose.graph[1][end]["O"] in pose.graph[1][end]["C"].children
        @test pose.graph[1][end]["OXT"] === nothing

        ProtoSyn.Peptides.uncap!(pose, rid"1")
        @test pose.state.size == 38
        @test length(pose.graph[1][1]["N"].children) == 1
        @test length(pose.graph[1][end]["C"].children) == 1
        @test pose.graph[1][1]["H"] === nothing
        @test pose.graph[1][1]["H1"] === nothing
        @test pose.graph[1][1]["H2"] === nothing
        @test pose.graph[1][1]["H3"] === nothing
        @test pose.graph[1][end]["O"] in pose.graph[1][end]["C"].children
        @test pose.graph[1][end]["OXT"] === nothing
    end

    @testset verbose = true "$(@sprintf "%-54s" "Cap (no selection)")" begin
        pose = copy(backup)
        @test pose.state.size == 39
        @test length(pose.graph[1][1]["N"].children) == 2
        @test length(pose.graph[1][end]["C"].children) == 1
        @test pose.graph[1][1]["H"] in pose.graph[1][1]["N"].children
        @test pose.graph[1][1]["H1"] === nothing
        @test pose.graph[1][1]["H2"] === nothing
        @test pose.graph[1][1]["H3"] === nothing
        @test pose.graph[1][end]["O"] in pose.graph[1][end]["C"].children
        @test pose.graph[1][end]["OXT"] === nothing

        ProtoSyn.Peptides.cap!(pose)
        @test pose.state.i2c == false
        @test pose.state.c2i == false
        @test pose.state.size == 42
        @test length(pose.graph[1][1]["N"].children) == 4
        @test length(pose.graph[1][end]["C"].children) == 2
        @test pose.graph[1][1]["H"] === nothing
        @test pose.graph[1][1]["H1"] in pose.graph[1][1]["N"].children
        # @test collect(pose.state[pose.graph[1][1]["H1"]].t) ≈ [0.22963328369701852, -1.7258132932668198, 0.697509161768422] atol = 1e-5
        @test pose.graph[1][1]["H2"] in pose.graph[1][1]["N"].children
        @test pose.graph[1][1]["H3"] in pose.graph[1][1]["N"].children
        @test pose.graph[1][end]["O"] in pose.graph[1][end]["C"].children
        @test pose.graph[1][end]["OXT"] in pose.graph[1][end]["C"].children
        @test collect(pose.state[pose.graph[1][end]["OXT"]].t) ≈ [7.637536790481703, -7.4116900395652, -0.4484369291709493] atol = 1e-5
    end

    @testset verbose = true "$(@sprintf "%-54s" "Cap (with selection)")" begin
        pose = copy(backup)
        @test pose.state.size == 39
        @test length(pose.graph[1][1]["N"].children) == 2
        @test length(pose.graph[1][end]["C"].children) == 1
        @test pose.graph[1][1]["H"] in pose.graph[1][1]["N"].children
        @test pose.graph[1][1]["H1"] === nothing
        @test pose.graph[1][1]["H2"] === nothing
        @test pose.graph[1][1]["H3"] === nothing
        @test pose.graph[1][end]["O"] in pose.graph[1][end]["C"].children
        @test pose.graph[1][end]["OXT"] === nothing

        ProtoSyn.Peptides.cap!(pose, rid"1")
        @test pose.state.i2c == false
        @test pose.state.c2i == false
        @test pose.state.size == 41
        @test length(pose.graph[1][1]["N"].children) == 4
        @test length(pose.graph[1][end]["C"].children) == 1
        @test pose.graph[1][1]["H"] === nothing
        @test pose.graph[1][1]["H1"] in pose.graph[1][1]["N"].children
        # @test collect(pose.state[pose.graph[1][1]["H1"]].t) ≈ [0.22963329021971224, -1.7472235478195153, 0.6757668454837074] atol = 1e-5
        @test pose.graph[1][1]["H2"] in pose.graph[1][1]["N"].children
        @test pose.graph[1][1]["H3"] in pose.graph[1][1]["N"].children
        @test pose.graph[1][end]["O"] in pose.graph[1][end]["C"].children
        @test pose.graph[1][end]["OXT"] === nothing
    end
end