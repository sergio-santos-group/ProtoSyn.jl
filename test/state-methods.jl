@testset "Measures" begin
    @test ProtoSyn.distance(pose.state[1], pose.state[2]) ≈ 1.0093879999999997 atol = 1e-10
    @test ProtoSyn.angle(pose.state[1], pose.state[2], pose.state[3]) ≈ 0.4174336256191947 atol = 1e-10
    @test ProtoSyn.dihedral(pose.state[1], pose.state[2], pose.state[3], pose.state[4]) ≈ -1.1582374633384385 atol = 1e-10
end