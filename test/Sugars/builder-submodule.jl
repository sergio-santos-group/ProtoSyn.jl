@testset verbose = true "Materials methods" begin

    @testset verbose = true "Grammar" begin
        g = ProtoSyn.Sugars.grammar("amylose");
        @test isa(g, ProtoSyn.LGrammar)
        @test length(g.operators) === 4
        @test length(g.variables) === 3
    end

    @testset verbose = true "Builder" begin
        g = ProtoSyn.Sugars.grammar("amylose");
        pose = ProtoSyn.build(g, seq"AAAβB[ɣCɣCɣC]AAA")
        @test pose.state.size === 206
        @test length(collect(eachatom(pose.graph))) === 206
        @test length(ProtoSyn.root(pose.graph).children) === 1
        @test length(ProtoSyn.root(pose.graph).container.children) === 1
        @test length(pose.graph[1][3].children) === 1
        @test length(pose.graph[1][4].children) === 2
        @test pose.graph[1][4]["O6"].parent == pose.graph[1][3]["C1"]
        @test pose.graph[1][5]["C1"] in pose.graph[1][4]["O4"].children
        @test pose.graph[1][8]["O4"] in pose.graph[1][4]["C1"].children
        @test length(pose.graph[1][5].children) === 1
    end
end