@testset verbose = true "Materials methods" begin

    @testset verbose = true "Primitive" begin
        pose = ProtoSyn.Materials.primitive()
        @test pose.state.size === 1
        @test length(pose.state.items) === 4
        @test length(collect(eachatom(pose.graph))) === 1
        @test pose.graph[1][1][1].name == "C"
    end

    @testset verbose = true "Body centered" begin
        pose = ProtoSyn.Materials.body_centered()
        @test pose.state.size === 2
        @test length(pose.state.items) === 5
        @test length(collect(eachatom(pose.graph))) === 2
        @test pose.graph[1][1][1].name == "C1"
        @test pose.graph[1][1][2].name == "C2"
    end
end