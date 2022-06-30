using StaticArrays

@testset verbose = true "Mutators $(repeat("-", 48))" begin

    @testset verbose = true "$(@sprintf "%-54s" "Dihedral")" begin
        pose = copy(backup)
        m = ProtoSyn.Mutators.DihedralMutator(() -> pi, 1.0, 1.0, an"C" & rid"2")
        ∠1_1 = ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["C"])
        ∠1_2 = ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CA"])
        ∠1_3 = ProtoSyn.getdihedral(pose.state, pose.graph[1][3]["N"])
        m(pose); sync!(pose)
        ∠2_1 = ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["C"])
        ∠2_2 = ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CA"])
        ∠2_3 = ProtoSyn.getdihedral(pose.state, pose.graph[1][3]["N"])
        @test ∠2_1  ≈ ∠1_1 + π atol = 1e-5
        @test ∠2_2  === ∠1_2
        @test ∠2_3  === ∠1_3
    end

    @testset verbose = true "$(@sprintf "%-54s" "Crankshaft")" begin
        pose = copy(backup)
        m = ProtoSyn.Mutators.CrankshaftMutator(()->deg2rad(5), 1.0, 1.0, an"CA" & (rid"1" | rid"2"), !(an"^CA$|^N$|^C$|^H$|^O$"r))
        ∠1_1 = ProtoSyn.getdihedral(pose.state, pose.graph[1][1]["C"])
        ∠1_2 = ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CA"])
        ∠1_3 = ProtoSyn.getdihedral(pose.state, pose.graph[1][3]["N"])
        CB_pos_1 = pose.state[pose.graph[1][2]["CB"]].t
        m(pose); sync!(pose)
        ∠2_1 = ProtoSyn.getdihedral(pose.state, pose.graph[1][1]["C"])
        ∠2_2 = ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CA"])
        ∠2_3 = ProtoSyn.getdihedral(pose.state, pose.graph[1][3]["N"])
        CB_pos_2 = pose.state[pose.graph[1][2]["CB"]].t
        @test ∠2_1  !== ∠1_1
        @test ∠2_1  ≈ -3.1156004909929313 atol = 1e-5
        @test ∠2_3  !== ∠1_3
        @test ∠2_3  ≈ 3.1154832277723146 atol = 1e-5
        @test ∠1_2 === ∠2_2
        @test CB_pos_1 !== CB_pos_2
        T = eltype(pose.state)
        @test CB_pos_2 ≈ [3.8810175757721384, -4.721163091936314, -1.3002006003050992] atol = 1e-5
    end

    @testset verbose = true "$(@sprintf "%-54s" "Translation Rigid Body")" begin
        pose = copy(backup)
        m = ProtoSyn.Mutators.TranslationRigidBodyMutator(()->[1.0, 1.0, 1.0], 1.0, nothing)
        pos_1_1 = pose.state[pose.graph[1][1]["CA"]].t
        pos_1_2 = pose.state[pose.graph[1][2]["CA"]].t
        pos_1_3 = pose.state[pose.graph[1][3]["CA"]].t
        m(pose); sync!(pose)
        pos_2_1 = pose.state[pose.graph[1][1]["CA"]].t
        pos_2_2 = pose.state[pose.graph[1][2]["CA"]].t
        pos_2_3 = pose.state[pose.graph[1][3]["CA"]].t
        @test pos_1_1 .+ [1.0, 1.0, 1.0] == pos_2_1
        @test pos_1_2 .+ [1.0, 1.0, 1.0] == pos_2_2
        @test pos_1_3 .+ [1.0, 1.0, 1.0] == pos_2_3
    end

    @testset verbose = true "$(@sprintf "%-54s" "Rotation Rigid Body")" begin
        pose = copy(backup)
        m = ProtoSyn.Mutators.RotationRigidBodyMutator(()->[1.0, 1.0, 1.0], ()->deg2rad(10), ProtoSyn.center_of_mass, 1.0, nothing)
        pos_1_1 = pose.state[pose.graph[1][1]["CA"]].t
        pos_1_2 = pose.state[pose.graph[1][2]["CA"]].t
        pos_1_3 = pose.state[pose.graph[1][3]["CA"]].t
        m(pose); sync!(pose)
        pos_2_1 = pose.state[pose.graph[1][1]["CA"]].t
        pos_2_2 = pose.state[pose.graph[1][2]["CA"]].t
        pos_2_3 = pose.state[pose.graph[1][3]["CA"]].t
        @test pos_1_1 !== pos_2_1
        @test pos_1_2 !== pos_2_2
        @test pos_1_3 !== pos_2_3
        T = eltype(pose.state)
        @test pos_2_1 ≈ [1.7573120826160027, -1.4568106376562309, 0.7092469773463803]
        @test pos_2_2 ≈ [4.237565726518101, -4.178262545553745, 0.18364182332261403]
        @test pos_2_3 ≈ [7.863391631896685, -4.874122906544416, -0.2690548282792112]
    end

    @testset verbose = true "$(@sprintf "%-54s" "Backrub")" begin
        pose = copy(backup)
        m =  ProtoSyn.Mutators.BackrubMutator(()->[1.0, 1.0, 1.0], 1.0, 1.0, aid"2")
        pos_1_2 = pose.state[pose.graph[1][1][2]].t
        m(pose); sync!(pose)
        @test pose.state[pose.graph[1][1][2]].t == pos_1_2 .+ 1.0
    end

    @testset verbose = true "$(@sprintf "%-54s" "Compound")" begin
        pose = copy(backup)
        m1 = ProtoSyn.Mutators.DihedralMutator(() -> pi, 1.0, 1.0, an"C")
        m2 = ProtoSyn.Mutators.DihedralMutator(() -> pi, 1.0, 1.0, an"C")
        m =  ProtoSyn.Mutators.CompoundMutator([m1, m2], rid"2")
        ∠1_1 = ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["C"])
        ∠1_2 = ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CA"])
        ∠1_3 = ProtoSyn.getdihedral(pose.state, pose.graph[1][3]["N"])
        m(pose); sync!(pose)
        ∠2_1 = ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["C"])
        ∠2_2 = ProtoSyn.getdihedral(pose.state, pose.graph[1][2]["CA"])
        ∠2_3 = ProtoSyn.getdihedral(pose.state, pose.graph[1][3]["N"])
        @test ProtoSyn.unit_circle(∠2_1)  ≈ ProtoSyn.unit_circle(∠1_1 + 2π) atol = 1e-5
        @test ProtoSyn.unit_circle(∠2_2)  ≈ ProtoSyn.unit_circle(∠1_2 + 2π) atol = 1e-5
        @test ProtoSyn.unit_circle(∠2_3)  ≈ ProtoSyn.unit_circle(∠1_3 + 2π) atol = 1e-5
    end
end