@testset verbose = true "Pose methods $(repeat("-", 44))" begin

    @testset verbose = true "Aligning two poses" begin
        pose1 = copy(backup)
        pose2 = copy(backup)

        @test ProtoSyn.align!(pose, pose2) === pose
        @test ProtoSyn.rmsd(pose, pose2) ≈ 0.0 atol = 1e-10
    end

    @testset verbose = true "Pop atom (H)" begin
        pose = copy(backup)

        @test pose.state.size == 39
        @test pose.graph[1][2].size == 17
        @test length(pose.graph[1][2].items) == 17
        @test pose.graph[1][2]["H"] in pose.graph[1][2]["N"].children
        @test pose.graph[1][2]["H"] in pose.graph[1][2]["N"].bonds
        @test pose.graph[1][2].parent == pose.graph[1][1]

        ProtoSyn.pop_atom!(pose, pose.graph[1][2]["H"])
        @test pose.graph[1][2]["H"] === nothing
        @test pose.state.size == 38
        @test pose.graph[1][2].size == 16
        @test length(pose.graph[1][2].items) == 16
        @test !(pose.graph[1][2]["H"] in pose.graph[1][2]["N"].children)
        @test !(pose.graph[1][2]["H"] in pose.graph[1][2]["N"].bonds)
        @test pose.graph[1][2].parent == pose.graph[1][1]
        @test pose.state.i2c == false
    end

    @testset verbose = true "Pop atom (N)" begin
        pose = copy(backup)

        @test pose.state.size == 39
        @test pose.graph[1][2].size == 17
        @test length(pose.graph[1][2].items) == 17
        @test pose.graph[1][2]["N"] in pose.graph[1][1]["C"].children
        @test pose.graph[1][2]["N"] in pose.graph[1][1]["C"].bonds
        @test pose.graph[1][2]["N"] in pose.graph[1][2]["CA"].bonds
        @test pose.graph[1][2]["CA"].parent == pose.graph[1][2]["N"]
        @test pose.graph[1][2].parent == pose.graph[1][1]
        @test pose.graph[1][2]["CA"].id == 10
        @test pose.graph[1][2]["CA"].index == 10
        @test pose.graph[1][2]["CA"].ascendents == (10, 8, 6, 3)

        ProtoSyn.pop_atom!(pose, pose.graph[1][2]["N"])
        @test pose.graph[1][2]["N"] === nothing
        @test pose.state.size == 38
        @test pose.graph[1][2].size == 16
        @test length(pose.graph[1][2].items) == 16
        @test !(pose.graph[1][2]["N"] in pose.graph[1][1]["C"].children)
        @test !(pose.graph[1][2]["N"] in pose.graph[1][1]["C"].bonds)
        @test !(pose.graph[1][2]["N"] in pose.graph[1][2]["CA"].bonds)
        @test pose.graph[1][2]["CA"].parent == ProtoSyn.root(pose.graph)
        @test pose.graph[1][2].parent == ProtoSyn.root(pose.graph).container
        @test pose.graph[1][2]["H"].parent == ProtoSyn.root(pose.graph)
        @test pose.state.i2c == false
        @test pose.graph[1][2]["CA"].id == 9
        @test pose.graph[1][2]["CA"].index == 9
        @test pose.graph[1][2]["CA"].ascendents == (9, 0, -1, -2)
    end

    @testset verbose = true "Pop atom (first atom of pose)" begin
        pose = copy(backup)

        @test pose.state.size == 39
        @test pose.graph[1][1].size == 7
        @test length(pose.graph[1][1].items) == 7
        @test pose.graph[1][1]["N"] in ProtoSyn.root(pose.graph).children
        @test pose.graph[1][1]["N"] in pose.graph[1][1]["CA"].bonds
        @test pose.graph[1][1]["CA"].parent == pose.graph[1][1]["N"]
        @test pose.graph[1][1].parent == ProtoSyn.root(pose.graph).container
        @test pose.graph[1][1]["CA"].id == 3
        @test pose.graph[1][1]["CA"].index == 3
        @test pose.graph[1][1]["CA"].ascendents == (3, 1, 0, -1)

        ProtoSyn.pop_atom!(pose, pose.graph[1][1]["N"])
        @test pose.graph[1][1]["N"] === nothing
        @test pose.state.size == 38
        @test pose.graph[1][1].size == 6
        @test length(pose.graph[1][1].items) == 6
        @test !(pose.graph[1][1]["N"] in ProtoSyn.root(pose.graph).children)
        @test !(pose.graph[1][1]["N"] in pose.graph[1][1]["CA"].bonds)
        @test pose.graph[1][1]["CA"].parent == ProtoSyn.root(pose.graph)
        @test pose.graph[1][1].parent == ProtoSyn.root(pose.graph).container
        @test pose.graph[1][1]["H"].parent == ProtoSyn.root(pose.graph)
        @test pose.state.i2c == false
        @test pose.graph[1][1]["CA"].id == 2
        @test pose.graph[1][1]["CA"].index == 2
        @test pose.graph[1][1]["CA"].ascendents == (2, 0, -1, -2)
    end

    @testset verbose = true "Remove hydrogens" begin
        pose = copy(backup)

        @test pose.state.size == 39
        @test pose.graph[1][1].size == 7
        @test length(pose.graph[1][1].items) == 7
        @test pose.graph[1][1]["N"] in ProtoSyn.root(pose.graph).children
        @test pose.graph[1][1]["N"] in pose.graph[1][1]["CA"].bonds
        @test pose.graph[1][1]["CA"].parent == pose.graph[1][1]["N"]
        @test pose.graph[1][1].parent == ProtoSyn.root(pose.graph).container
        @test pose.graph[1][1]["CA"].id == 3
        @test pose.graph[1][1]["CA"].index == 3
        @test pose.graph[1][1]["CA"].ascendents == (3, 1, 0, -1)

        ProtoSyn.pop_atoms!(pose, as"H")

        @test pose.state.size == 21
        @test count(as"H"(pose)) === 0
        @test pose.graph[1][1].size == 4
        @test length(pose.graph[1][1].items) == 4
        @test !("H" in [a.name for a in pose.graph[1][1]["N"].bonds])
        @test !("H" in [a.name for a in pose.graph[1][1]["N"].children])
        @test pose.state.i2c == false
        @test pose.graph[1][1]["CA"].id == 2
        @test pose.graph[1][1]["CA"].index == 2
        @test pose.graph[1][1]["CA"].ascendents == (2, 1, 0, -1)
    end

    @testset verbose = true "Remove hydrogens" begin
        pose = copy(backup)
        ProtoSyn.pop_atoms!(pose, as"H")

        @test pose.state.size == 21
        @test count(as"H"(pose)) === 0
        @test pose.graph[1][1].size == 4
        @test length(pose.graph[1][1].items) == 4
        @test !("H" in [a.name for a in pose.graph[1][1]["N"].bonds])
        @test !("H" in [a.name for a in pose.graph[1][1]["N"].children])
        @test pose.state.i2c == false
        @test pose.graph[1][1]["CA"].id == 2
        @test pose.graph[1][1]["CA"].index == 2
        @test pose.graph[1][1]["CA"].ascendents == (2, 1, 0, -1)

        ProtoSyn.add_hydrogens!(pose, ProtoSyn.Peptides.grammar, nothing)

        @test pose.state.size == 39
        @test pose.graph[1][1].size == 7
        @test length(pose.graph[1][1].items) == 7
        @test pose.graph[1][1]["N"] in ProtoSyn.root(pose.graph).children
        @test pose.graph[1][1]["N"] in pose.graph[1][1]["CA"].bonds
        @test pose.graph[1][1]["CA"].parent == pose.graph[1][1]["N"]
        @test pose.graph[1][1].parent == ProtoSyn.root(pose.graph).container
        @test pose.graph[1][1]["CA"].id == 3
        @test pose.graph[1][1]["CA"].index == 3
        @test pose.graph[1][1]["CA"].ascendents == (3, 1, 0, -1)
    end
end