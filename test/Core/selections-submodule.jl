println("-----------\n Selections submodule:")

@testset "SerialSelection" begin
    pose = copy(backup)
    
    @test count(aid"1"(pose)) === 1
    @test count(rid"1"(pose)) === 1
end

@testset "RangeSelection" begin
    pose = copy(backup)
    
    @test count(aid"1:10"(pose)) === 10
    @test count(rid"1:3"(pose)) === 3
    @test count(rid"1:10"(pose)) === 3
end

@testset "FieldSelection" begin
    pose = copy(backup)
    
    @test count(an"CA"(pose)) === 3
    @test count(rn"GLY"(pose)) === 1
    @test count(rn"GLN"(pose)) === 0
end

@testset "TerminalSelection" begin
    pose = copy(backup)
    
    @test TerminalSelection()(pose).content == [true, false, true]
end

@testset "DistanceSelection" begin
    pose = copy(backup)
    
    @test count((5.0:aid"1")(pose)) === 12
    @test count((5.0:rid"1")(pose)) === 22
end

@testset "RandomSelection" begin
    pose = copy(backup)
    
    @test count(RandomSelection{Atom}()(pose)) === 1
    @test count(RandomSelection{Residue}()(pose)) === 1
end

@testset "RandomRangeSelection" begin
    pose = copy(backup)
    
    @test count(RandomRangeSelection{Atom}()(pose)) > 0
end

@testset "TrueSelection" begin
    pose = copy(backup)
    
    @test count(TrueSelection{Atom}()(pose)) == size(pose.graph)[end]
end

@testset "UnarySelection" begin
    pose = copy(backup)
    
    @test count(!an"CA"(pose)) == size(pose.graph)[end] - size(pose.graph)[end - 1]
end

@testset "Promotion" begin
    pose = copy(backup)
    
    @test ProtoSyn.selection_type(ProtoSyn.promote(an"CA", Residue)) === Residue
    @test typeof(ProtoSyn.gather(ProtoSyn.promote(an"CA"(pose), Residue, pose.graph), pose.graph)[1]) === Residue
end

@testset "BinarySelection" begin
    pose = copy(backup)
    
    @test count((an"CA" & an"CB")(pose)) === 0
    @test count((an"CA" | an"CB")(pose)) === 5
    @test count((an"CA" & (an"CA" | an"CB"))(pose)) === 3
    @test count((an"CA" & an"CA" | an"CB")(pose)) === 5
end