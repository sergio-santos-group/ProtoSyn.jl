@testset verbose = true "Peptides | Selections $(repeat("-", 35))" begin
    @testset verbose = true "$(@sprintf "%-54s" "Polar Selection")" begin
        pose = copy(backup)
        @test PolarSelection()(pose).content == [false, false, true]
    end

    @testset verbose = true "$(@sprintf "%-54s" "Sidechain Selection")" begin
        pose = copy(backup)
        @test count(SidechainSelection()(pose)) == 24
    end

    @testset verbose = true "$(@sprintf "%-54s" "Secondary Structure Selection")" begin
        pose = copy(backup)
        ProtoSyn.Peptides.setss!(pose, ProtoSyn.Peptides.SecondaryStructure[:helix])
        @test ss"helix"(pose).content == [true, true, false]
        @test SecondaryStructureSelection(:helix, deg2rad(50))(pose).content == [true, true, false]
    end

    @testset verbose = true "$(@sprintf "%-54s" "Protein Selection")" begin
        pose = copy(backup)
        @test ProteinSelection()(pose).content == [true, true, true]
    end

    @testset verbose = true "$(@sprintf "%-54s" "Chi Selection")" begin
        pose = copy(backup)
        @test count(chi"1"(pose.graph[1, 1])) === 0
        @test count(chi"1"(pose.graph[1, 2])) === 1
        @test count(chi"2"(pose.graph[1, 2])) === 1
        @test count(chi"3"(pose.graph[1, 2])) === 1
        @test count(chi"4"(pose.graph[1, 2])) === 0

        @test count(ChiSelection(1)(pose.graph[1, 1])) === 0
        @test count(ChiSelection(1)(pose.graph[1, 2])) === 1
        @test count(ChiSelection(2)(pose.graph[1, 2])) === 1
        @test count(ChiSelection(3)(pose.graph[1, 2])) === 1
        @test count(ChiSelection(4)(pose.graph[1, 2])) === 0
    end

    @testset verbose = true "$(@sprintf "%-54s" "Phi/Psi/Omega Selection")" begin
        pose = copy(backup)

        @test count(PhiSelection()(pose.graph[1, 1])) === 0
        @test count(PhiSelection()(pose.graph[1, 2])) === 1
        @test count(PhiSelection()(pose.graph[1, 3])) === 1
        @test count(PsiSelection()(pose.graph[1, 1])) === 1
        @test count(PsiSelection()(pose.graph[1, 2])) === 1
        @test count(PsiSelection()(pose.graph[1, 3])) === 0
        @test count(OmegaSelection()(pose.graph[1, 1])) === 0
        @test count(OmegaSelection()(pose.graph[1, 2])) === 1
        @test count(OmegaSelection()(pose.graph[1, 3])) === 1
    end
end