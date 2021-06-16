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
end