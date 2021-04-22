using ProtoSyn.Peptides

println("-----------\n Builder:")

res_lib = Peptides.grammar(Float64)
pose    = ProtoSyn.build(res_lib, seq"GME")

@testset "Copying a pose" begin
    backup = copy(pose)
    @test size(pose.graph) === size(backup.graph)
end

backup = copy(pose)

@testset "Appending from derivation" begin
    pose = copy(backup)
    ProtoSyn.append_fragment!(pose, pose.graph[1][3], res_lib, seq"MMM")
    @test pose.graph[1][4] in pose.graph[1][3].children
    @test pose.graph[1][4]["N"] in pose.graph[1][3]["C"].children
    @test ProtoSyn.count_residues(pose.graph) === 6
end

@testset "Unbonding two atoms" begin
    pose = copy(backup)
    ProtoSyn.unbond(pose, pose.graph[1][1]["C"], pose.graph[1, 2, "N"])
    @test !(pose.graph[1][2]["N"] in pose.graph[1][1]["C"].bonds)
    @test pose.graph[1][2].parent === root(pose.graph).container
    @test pose.graph[1][2]["N"].parent === root(pose.graph)
end

@testset "Inserting from derivation" begin
    pose = copy(backup)
    ProtoSyn.unbond(pose, pose.graph[1][1]["C"], pose.graph[1, 2, "N"])
    ProtoSyn.insert_fragment!(pose, pose.graph[1][2], res_lib, seq"MMM")
    @test pose.graph[1][2] in pose.graph[1][1].children
    @test pose.graph[1][2]["N"] in pose.graph[1][1]["C"].children
    @test ProtoSyn.count_residues(pose.graph) === 6
end