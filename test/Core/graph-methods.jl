@testset verbose = true "Graph methods $(repeat("-", 43))" begin

    @testset verbose = true "Root & Origin" begin
        pose = copy(backup)

        @test ProtoSyn.root(pose.graph) === pose.graph[1][1][1].parent
        @test ProtoSyn.origin(pose.graph) === pose.graph[1][1][1]
    end

    @testset verbose = true "Parenthood relationships" begin
        pose = copy(backup)

        @test hasparent(pose.graph[1][1][1])
        @test hasparent(pose.graph[1][1])
        @test isparent(pose.graph[1][1], pose.graph[1][2])
        @test isparent(ProtoSyn.root(pose.graph), pose.graph[1][1][1])
        @test isparent(ProtoSyn.root(pose.graph).container, pose.graph[1][1])
        @test_throws ErrorException setparent!(pose.graph[1][1][1], pose.graph[1][2][1])
        @test_throws ErrorException  setparent!(pose.graph[1][1], pose.graph[1][2])
        @test popparent!(pose.graph[1][2][1]) === pose.graph[1][2][1]
        @test !hasparent(pose.graph[1][2][1])
        @test setparent!(pose.graph[1][2][1], pose.graph[1][1][1]) === pose.graph[1][2][1]
        @test hasparent(pose.graph[1][2][1])
        @test ascendents(pose.graph[1][2][1], 4) === (8, 1, 0, -1)
    end

    @testset verbose = true "Detach" begin
        pose = copy(backup)

        s = detach(pose.graph[1])
        @test size(pose.graph) === (0, 0, 0)
        @test typeof(s) === Segment
        @test origin(s) === s[1][1]
        @test root(s) === nothing
        @test length(root(pose.graph).children) === 0
    end

    @testset verbose = true "Graph travel" begin
        pose = copy(backup)

        @test length(travel_graph(pose.graph[1][3][1])) === 15
        @test length(unique([a.container.id for a in travel_graph(pose.graph[1][3][1])])) === 1
        @test ProtoSyn.travel_graph(pose.graph[1, 3, "CA"], search_algorithm = ProtoSyn.BFS)[1:4] == [pose.graph[1, 3, "CA"], pose.graph[1, 3, "HA"], pose.graph[1, 3, "C"], pose.graph[1, 3, "CB"]]
        @test ProtoSyn.travel_graph(pose.graph[1, 3, "CA"], search_algorithm = ProtoSyn.DFS)[1:4] == [pose.graph[1, 3, "CA"], pose.graph[1, 3, "CB"], pose.graph[1, 3, "CG"], pose.graph[1, 3, "CD"]]
        @test ProtoSyn.is_contiguous(pose, rid"2:3")
        @test ProtoSyn.is_contiguous(pose, an"K") === nothing
        @test !ProtoSyn.is_contiguous(pose, rid"1" | rid"3")
    end

    @testset verbose = true "Container manipulation" begin
        pose = copy(backup)

        @test hascontainer(pose.graph[1][1][1])
        @test hascontainer(pose.graph[1][1])
        @test !hascontainer(pose.graph)
        @test push!(pose.graph[1][2], pose.graph[1][1][1]) === pose.graph[1][2]
        @test length(pose.graph[1][2]) === 18
        @test insert!(pose.graph[1][2], 1, pose.graph[1][1][1]) === pose.graph[1][2]
        @test length(pose.graph[1][2]) === 19
        @test delete!(pose.graph[1][2], pose.graph[1][2][2]) === pose.graph[1][2]
        @test length(pose.graph[1][2]) === 18
    end

    @testset verbose = true "Indexation" begin
        pose = copy(backup)
        pose.graph[1][1][1].id = 1000

        @test typeof(genid()) === Int64
        @test reindex(pose.graph) === pose.graph
        @test pose.graph[1][1][1].id === 1
        @test ids(pose.graph[1][1].items) ==  [1, 2, 3, 4, 5, 6, 7]
    end

    @testset verbose = true "Unbonding (intra-residue)" begin
        pose = copy(backup)

        @test unbond!(pose, pose.graph[1][3]["N"], pose.graph[1][3]["H"]) === pose
        @test !(pose.graph[1][3]["H"] in pose.graph[1][3]["N"].bonds)
        @test pose.graph[1][3]["H"].parent == ProtoSyn.root(pose.graph)
        @test pose.graph[1][3].parent == pose.graph[1][2]
        @test pose.state.i2c == true
    end

    @testset verbose = true "Unbonding (inter-residues)" begin
        pose = copy(backup)

        @test unbond!(pose, pose.graph[1][2]["C"], pose.graph[1][3]["N"]) === pose
        @test !(pose.graph[1][2]["C"] in pose.graph[1][3]["N"].bonds)
        @test pose.graph[1][3]["N"].parent == ProtoSyn.root(pose.graph)
        @test pose.graph[1][3].parent == ProtoSyn.root(pose.graph).container
        @test pose.state.i2c == true
    end

    @testset verbose = true "Re-bonding" begin
        pose = copy(backup)

        @test unbond!(pose, pose.graph[1][1]["C"], pose.graph[1][2]["N"]) === pose
        @test !(pose.graph[1][2]["N"] in pose.graph[1][1]["C"].bonds)
        ProtoSyn.bond(pose.graph[1][2]["N"], pose.graph[1][1]["C"])
        @test pose.graph[1][2]["N"] in pose.graph[1][1]["C"].bonds
        @test ProtoSyn.bond(pose.graph[1][2]["N"], pose.graph[1][1]["C"]) === nothing
    end

    @testset verbose = true "Unbonding + Dihedral rotation" begin
        pose = copy(backup)
        
        @test length(pose.graph[1][2]["CB"].children) == 3
        @test pose.graph[1][2]["CG"].parent == pose.graph[1][2]["CB"]
        ProtoSyn.setdihedral!(pose.state, pose.graph[1][2]["CG"], deg2rad(90))
        ProtoSyn.unbond!(pose, pose.graph[1][2]["CG"], pose.graph[1][2]["CB"])
        @test pose.graph[1][2]["CG"].parent == ProtoSyn.root(pose.graph)
        @test pose.state.i2c == true
        sync!(pose)
        @test pose.state.i2c == false
        @test length(pose.graph[1][2]["CB"].children) == 2
        @test all(collect(pose.state[pose.graph[1][2]["CG"]].t) .≈ [4.800701090587441, -4.607782270435299, -2.3619776332006834])
    end

end