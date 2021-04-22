@testset "Aligning two poses" begin
    pose1 = copy(backup)
    pose2 = copy(backup)

    @test ProtoSyn.align!(pose, pose2) === pose
    @test ProtoSyn.rmsd(pose, pose2) ≈ 0.0 atol = 1e-10
end