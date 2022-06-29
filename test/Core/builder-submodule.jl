using ProtoSyn.Peptides

@testset verbose = true "Core | Builder submodule $(repeat("-", 32))" begin

    @testset verbose = true "Building a pose from derivation" begin
        res_lib = Peptides.grammar
        @test length(res_lib.variables) === 21
        @test length(res_lib.operators) === 3
        frag    = fragment(res_lib, seq"GME") 
        @test size(frag.graph) === (3, 39)
        @test typeof(frag.graph) === Segment
        pose    = ProtoSyn.build(res_lib, seq"GME")
        @test size(pose.graph) === (1, 3, 39)
        @test typeof(pose.graph) === Topology
        @test length(pose.state) === pose.state.size === 39
        @test pose.state[1].t == pose.state.x[:, 1]
    end

    res_lib = Peptides.grammar
    pose    = ProtoSyn.build(res_lib, seq"GME")
    sync!(pose)

    @testset verbose = true "Copying a pose" begin
        backup = copy(pose)
        @test size(pose.graph) === size(backup.graph)
    end

    global backup = copy(pose)

    @testset verbose = true "Appending from derivation" begin
        pose = copy(backup)
        ProtoSyn.append_fragment!(pose, pose.graph[1][3], res_lib, seq"MMM")
        @test pose.graph[1][4] in pose.graph[1][3].children
        @test pose.graph[1][4]["N"] in pose.graph[1][3]["C"].children
        @test ProtoSyn.count_residues(pose.graph) === 6
    end

    @testset verbose = true "Unbonding two atoms" begin
        pose = copy(backup)
        ProtoSyn.unbond!(pose, pose.graph[1][1]["C"], pose.graph[1, 2, "N"])
        @test !(pose.graph[1][2]["N"] in pose.graph[1][1]["C"].bonds)
        @test pose.graph[1][2].parent === ProtoSyn.root(pose.graph).container
        @test pose.graph[1][2]["N"].parent === ProtoSyn.root(pose.graph)
    end

    @testset verbose = true "Inserting from derivation" begin
        pose = copy(backup)
        ProtoSyn.unbond!(pose, pose.graph[1][1]["C"], pose.graph[1, 2, "N"])
        ProtoSyn.insert_fragment!(pose, pose.graph[1][2], res_lib, seq"MMM")
        @test !(pose.graph[1][2] in pose.graph[1][1].children)
        @test pose.graph[1][2].parent === ProtoSyn.root(pose.graph).container
        @test pose.graph[1][2]["N"].parent === ProtoSyn.root(pose.graph)
        @test pose.graph[1][3]["N"].parent === pose.graph[1][2]["C"]
        @test ProtoSyn.count_residues(pose.graph) === 6
    end

end