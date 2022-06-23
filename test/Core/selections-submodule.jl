@testset verbose = true "Selections submodule $(repeat("-", 36))" begin

    @testset verbose = true "SerialSelection     " begin
        pose = copy(backup)
        
        @test count(aid"1"(pose)) === 1
        @test count(rid"1"(pose)) === 1
    end

    @testset verbose = true "RangeSelection      " begin
        pose = copy(backup)
        
        @test count(aid"1:10"(pose)) === 10
        @test count(rid"1:3"(pose)) === 3
        @test count(rid"1:10"(pose)) === 3
    end

    @testset verbose = true "FieldSelection      " begin
        pose = copy(backup)
        
        @test count(an"CA"(pose)) === 3
        @test count(rn"GLY"(pose)) === 1
        @test count(rn"GLN"(pose)) === 0
    end

    @testset verbose = true "TerminalSelection   " begin
        pose = copy(backup)
        
        @test UpstreamTerminalSelection{Residue}()(pose).content == [true, false, false]
        @test DownstreamTerminalSelection{Residue}()(pose).content == [false, false, true]
        @test (DownstreamTerminalSelection{Residue}() | UpstreamTerminalSelection{Residue}())(pose).content == [true, false, true]
    end

    @testset verbose = true "DistanceSelection   " begin
        pose = copy(backup)
        
        @test count((5.0:aid"1")(pose)) === 12
        @test count((5.0:rid"1")(pose)) === 22
    end

    @testset verbose = true "RandomSelection     " begin
        pose = copy(backup)
        
        @test count(RandomSelection{Atom}()(pose)) === 1
        @test count(RandomSelection{Residue}()(pose)) === 1
    end

    @testset verbose = true "RandomRangeSelection" begin
        pose = copy(backup)
        
        @test count(RandomRangeSelection{Atom}()(pose)) > 0
    end

    @testset verbose = true "TrueSelection       " begin
        pose = copy(backup)
        
        @test count(TrueSelection{Atom}()(pose)) == size(pose.graph)[end]
    end

    @testset verbose = true "UnarySelection      " begin
        pose = copy(backup)
        
        @test count(!an"CA"(pose)) == size(pose.graph)[end] - size(pose.graph)[end - 1]
    end

    @testset verbose = true "Promotion           " begin
        pose = copy(backup)
        
        @test ProtoSyn.selection_type(ProtoSyn.promote(an"CA", Residue)) === Residue
        @test typeof(ProtoSyn.gather(ProtoSyn.promote(an"CA"(pose), Residue, pose.graph), pose.graph)[1]) === Residue
    end

    @testset verbose = true "BinarySelection     " begin
        pose = copy(backup)
        
        @test count((an"CA" & an"CB")(pose)) === 0
        @test count((an"CA" | an"CB")(pose)) === 5
        @test count((an"CA" & (an"CA" | an"CB"))(pose)) === 3
        @test count((an"CA" & an"CA" | an"CB")(pose)) === 5
    end

    @testset verbose = true "AromaticSelection     " begin
        pose = copy(backup)
        ProtoSyn.Peptides.mutate!(pose, pose.graph[1, 2], Peptides.grammar, seq"W")
        
        @test count(AromaticSelection()(pose)) === 9
        @test count(AromaticSelection(5)(pose)) === 5
        @test count(AromaticSelection(4)(pose)) === 0
    end

    @testset verbose = true "BondCountSelection     " begin
        pose = copy(backup)
        
        @test count(BondCountSelection(2)(pose)) === 3
        @test count(BondCountSelection(3, <)(pose)) === 26
        @test count(BondCountSelection(4, >=)(pose)) === 8
    end

    @testset verbose = true "BondedToSelection     " begin
        pose = copy(backup)
        
        @test count(BondedToSelection(an"C")(pose)) === 8
        @test count(BondedToSelection(an"CA")(pose)) === 12
    end

    @testset verbose = true "ChargeSelection     " begin
        pose = copy(backup)
        
        @test count(ChargeSelection(0.0, >=)(pose)) === 23
        @test count(ChargeSelection(0.0, <)(pose)) === 16
    end
end