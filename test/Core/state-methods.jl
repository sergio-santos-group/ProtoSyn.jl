@testset verbose = true "State methods $(repeat("-", 43))" begin

    @testset verbose = true "Cartesian & Internal Coordinates" begin
        pose = copy(backup)

        @test ProtoSyn.request_i2c!(pose.state).i2c
        @test !sync!(pose).state.i2c
        @test ProtoSyn.request_c2i!(pose.state).c2i
        @test !sync!(pose).state.c2i

        @test getdihedral(pose.state, pose.graph[1][2][1]) ≈ 3.141592653589793 atol = 1e-10
        @test setdihedral!(pose.state, pose.graph[1][2][1], deg2rad(10.0)) === pose.state
        @test pose.state.i2c
        @test ProtoSyn.i2c!(pose.state, pose.graph) === pose.state
        @test !pose.state.i2c
        @test getdihedral(pose.state, pose.graph[1][2][1]) ≈ 0.17453292519943275 atol = 1e-10
        @test rotate_dihedral!(pose.state, pose.graph[1][2][1], deg2rad(10)) === pose.state
        @test sync!(pose) === pose
        @test getdihedral(pose.state, pose.graph[1][2][1]) ≈ 0.3490658503988655 atol = 1e-10

        pose.state[10].t = zeros(eltype(pose.state), 3)
        @test ProtoSyn.request_c2i!(pose.state).c2i
        @test pose.state.c2i
        @test ProtoSyn.c2i!(pose.state, pose.graph) === pose.state
        @test !pose.state.c2i
        @test ProtoSyn.getdihedral(pose.state, pose.graph[1][2][3]) ≈ -0.1505983276216401 atol = 1e-10

        pose = copy(backup)
        @test ProtoSyn.setoffset!(pose.state, pose.graph[1][1][3], deg2rad(10)) === pose.state
        @test pose.state.i2c
        @test rad2deg(getdihedral(pose.state, pose.graph[1][1][3])) ≈ 10 atol = 1e-5
    end

    @testset verbose = true "Indexation" begin
        pose = copy(backup)
        pose.state[1].index = 10

        @test reindex(pose.state) === pose.state
        @test pose.state[1].index === 1
    end

    @testset verbose = true "Measures" begin
        pose = copy(backup)

        @test ProtoSyn.distance(pose.state[1], pose.state[2]) ≈ 1.0093879999999997 atol = 1e-10
        @test ProtoSyn.angle(pose.state[1], pose.state[2], pose.state[3]) ≈ 0.4174336256191947 atol = 1e-10
        @test ProtoSyn.dihedral(pose.state[1], pose.state[2], pose.state[3], pose.state[4]) ≈ -1.1582374633384385 atol = 1e-10
    end

end