println("-----------\n Calculators:")

@testset "Verlet List" begin
    vl1 = ProtoSyn.Calculators.VerletList(pose.state.size)
    vl2 = ProtoSyn.Calculators.VerletList(pose.state.size)
    vl3 = ProtoSyn.Calculators.VerletList(pose.state.size)
    @test vl1.size === pose.state.size
    vl1.cutoff = 8.0
    vl2.cutoff = 8.0
    vl3.cutoff = 8.0
    vl1 = ProtoSyn.Calculators.update!(ProtoSyn.SISD_0, vl1, pose)
    vl2 = ProtoSyn.Calculators.update!(ProtoSyn.SIMD_1, vl2, pose)
    @test vl1.offset == vl2.offset
    @test vl1.list[vl1.offset[end]] == vl2.list[vl2.offset[end]]
    vl3 = ProtoSyn.Calculators.update!(vl3, pose)
    @test vl1.offset == vl3.offset
    @test vl1.list[vl1.offset[end]] == vl3.list[vl3.offset[end]]

    vl4 = ProtoSyn.Calculators.VerletList(pose.state.size); vl4.cutoff = 8.0
    ProtoSyn.Calculators.update!(ProtoSyn.SISD_0, vl4, pose, an"CA")
    @test vl4.offset[1] === 0
    @test vl4.offset[3] === 1
    @test vl4.list[vl4.offset[3]] === 10

    vl5 = ProtoSyn.Calculators.VerletList(pose.state.size); vl5.cutoff = 8.0
    ProtoSyn.Calculators.update!(ProtoSyn.SIMD_1, vl5, pose, an"CA")

    @test vl4.offset == vl5.offset
    @test vl4.list[vl4.offset[10]] == vl5.list[vl5.offset[10]]
end


@testset "SISD_0 Distance Matrix" begin
    pose = copy(backup)
    dm = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SISD_0, pose.state.x.coords)

    @test dm[1, 1] === 0.0
    @test size(dm) === (39, 39)
    @test dm[end, 1] === 0.0
    @test dm[1, end] ≈ 9.503650710901654 atol = 1e-5
    @test ProtoSyn.Calculators.distance_matrix(ProtoSyn.SISD_0, pose.state) == dm
    @test ProtoSyn.Calculators.distance_matrix(ProtoSyn.SISD_0, pose) == dm
    
    dm_sele = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SISD_0, pose, an"CA")
    @test dm_sele[1, 1] === 0.0
    @test size(dm_sele) === (3, 3)
    @test dm_sele[end, 1] === 0.0
    @test dm_sele[1, end] ≈ 7.065359511597963 atol = 1e-5

    dm_full = ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.SISD_0, pose.state.x.coords)
    @test dm_full[1, 1] === 0.0
    @test size(dm_full) === (39, 39)
    @test dm_full[end, 1] === dm_full[1, end]
    @test dm_full[1, end] ≈ 9.503650710901654 atol = 1e-5
    @test ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.SISD_0, pose.state) == dm_full
    @test ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.SISD_0, pose) == dm_full

    dm_full_sele = ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.SISD_0, pose, an"CA")
    @test dm_full_sele[1, 1] === 0.0
    @test size(dm_full_sele) === (3, 3)
    @test dm_full_sele[end, 1] === dm_full_sele[1, end]
    @test dm_full_sele[1, end] ≈ 7.065359511597963 atol = 1e-5

    vl = ProtoSyn.Calculators.VerletList(pose.state.size); vl.cutoff = 8.0
    ProtoSyn.Calculators.update!(ProtoSyn.SISD_0, vl, pose, an"CA")
    dm_vl1 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SISD_0, pose.state.x.coords, vl)
    dm_vl2 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SISD_0, pose.state, vl)
    dm_vl3 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SISD_0, pose, vl)
    @test dm_vl1 == dm_vl2
    @test dm_vl1 == dm_vl3
    @test dm_vl1[3, 10] ≈ 3.7194380493597423 atol = 1e-5
end

@testset "SIMD_1 Distance Matrix" begin
    pose = copy(backup)
    dm = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state.x.coords)

    @test dm[1, 1] === 0.0
    @test size(dm) === (39, 39)
    @test dm[end, 1] === 0.0
    @test dm[1, end] ≈ 9.503650710901654 atol = 1e-5
    @test ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state) == dm
    @test ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose) == dm
    
    dm_sele = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose, an"CA")
    @test dm_sele[1, 1] === 0.0
    @test size(dm_sele) === (3, 3)
    @test dm_sele[end, 1] === 0.0
    @test dm_sele[1, end] ≈ 7.065359511597963 atol = 1e-5

    dm_full = ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.SIMD_1, pose.state.x.coords)
    @test dm_full[1, 1] === 0.0
    @test size(dm_full) === (39, 39)
    @test dm_full[end, 1] === dm_full[1, end]
    @test dm_full[1, end] ≈ 9.503650710901654 atol = 1e-5
    @test ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.SIMD_1, pose.state) == dm_full
    @test ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.SIMD_1, pose) == dm_full

    dm_full_sele = ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.SIMD_1, pose, an"CA")
    @test dm_full_sele[1, 1] === 0.0
    @test size(dm_full_sele) === (3, 3)
    @test dm_full_sele[end, 1] === dm_full_sele[1, end]
    @test dm_full_sele[1, end] ≈ 7.065359511597963 atol = 1e-5

    vl = ProtoSyn.Calculators.VerletList(pose.state.size); vl.cutoff = 8.0
    ProtoSyn.Calculators.update!(ProtoSyn.SIMD_1, vl, pose, an"CA")
    dm_vl1 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state.x.coords, vl)
    dm_vl2 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state, vl)
    dm_vl3 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose, vl)
    @test dm_vl1 == dm_vl2
    @test dm_vl1 == dm_vl3
    @test dm_vl1[3, 10] ≈ 3.7194380493597423 atol = 1e-5
end

@testset "CUDA_2 Distance Matrix" begin
    pose = copy(backup)
    dm = collect(ProtoSyn.Calculators.distance_matrix(ProtoSyn.CUDA_2, pose.state.x.coords))

    @test dm[1, 1] === 0.0
    @test size(dm) === (39, 39)
    @test dm[end, 1] === dm[1, end]
    @test dm[1, end] ≈ 9.503650710901654 atol = 1e-5
    @test ProtoSyn.Calculators.distance_matrix(ProtoSyn.CUDA_2, pose.state) == dm
    @test ProtoSyn.Calculators.distance_matrix(ProtoSyn.CUDA_2, pose) == dm
    
    dm_sele = collect(ProtoSyn.Calculators.distance_matrix(ProtoSyn.CUDA_2, pose, an"CA"))
    @test dm_sele[1, 1] === 0.0
    @test size(dm_sele) === (3, 3)
    @test dm_sele[end, 1] === dm_sele[1, end]
    @test dm_sele[1, end] ≈ 7.065359511597963 atol = 1e-5

    dm_full = collect(ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.CUDA_2, pose.state.x.coords))
    @test dm_full[1, 1] === 0.0
    @test size(dm_full) === (39, 39)
    @test dm_full[end, 1] === dm_full[1, end]
    @test dm_full[1, end] ≈ 9.503650710901654 atol = 1e-5
    @test ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.CUDA_2, pose.state) == dm_full
    @test ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.CUDA_2, pose) == dm_full

    dm_full_sele = collect(ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.CUDA_2, pose, an"CA"))
    @test dm_full_sele[1, 1] === 0.0
    @test size(dm_full_sele) === (3, 3)
    @test dm_full_sele[end, 1] === dm_full_sele[1, end]
    @test dm_full_sele[1, end] ≈ 7.065359511597963 atol = 1e-5

    vl = ProtoSyn.Calculators.VerletList(pose.state.size); vl.cutoff = 8.0
    ProtoSyn.Calculators.update!(ProtoSyn.SIMD_1, vl, pose, an"CA")
    dm_vl1 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state.x.coords, vl)
    dm_vl2 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state, vl)
    dm_vl3 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose, vl)
    @test dm_vl1 == dm_vl2
    @test dm_vl1 == dm_vl3
    @test dm_vl1[3, 10] ≈ 3.7194380493597423 atol = 1e-5
end

@testset "Dynamic Acceleration Type Distance Matrix" begin
    pose = copy(backup)
    dm = collect(ProtoSyn.Calculators.distance_matrix(pose.state.x.coords))

    @test dm[1, 1] === 0.0
    @test size(dm) === (39, 39)
    @test dm[1, end] ≈ 9.503650710901654 atol = 1e-5
    @test ProtoSyn.Calculators.distance_matrix(pose.state) == dm
    @test ProtoSyn.Calculators.distance_matrix(pose) == dm
    
    dm_sele = collect(ProtoSyn.Calculators.distance_matrix(pose, an"CA"))
    @test dm_sele[1, 1] === 0.0
    @test size(dm_sele) === (3, 3)
    @test dm_sele[1, end] ≈ 7.065359511597963 atol = 1e-5

    dm_full = collect(ProtoSyn.Calculators.full_distance_matrix(pose.state.x.coords))
    @test dm_full[1, 1] === 0.0
    @test size(dm_full) === (39, 39)
    @test dm_full[1, end] ≈ 9.503650710901654 atol = 1e-5
    @test ProtoSyn.Calculators.full_distance_matrix(pose.state) == dm_full
    @test ProtoSyn.Calculators.full_distance_matrix(pose) == dm_full

    dm_full_sele = collect(ProtoSyn.Calculators.full_distance_matrix(pose, an"CA"))
    @test dm_full_sele[1, 1] === 0.0
    @test size(dm_full_sele) === (3, 3)
    @test dm_full_sele[1, end] ≈ 7.065359511597963 atol = 1e-5
end

@testset "Potential" begin
    p = ProtoSyn.Calculators.get_flat_bottom_potential()
    @test p(100.0) === 0.0
    p = ProtoSyn.Calculators.get_flat_bottom_potential(d3=10.0, d4 = 20.0)
    @test p(100.0) === 1700.0

    @test BitArray([false true true;true false true;true true false]) == ProtoSyn.Calculators.intra_residue_mask(pose, an"CA").content
    girm = ProtoSyn.Calculators.get_intra_residue_mask(an"CA")
    @test BitArray([false true true;true false true;true true false]) == girm(pose).content

    @test BitArray([false true true;true false true;true true false]) == ProtoSyn.Calculators.diagonal_mask(pose, an"CA").content
    gdm = ProtoSyn.Calculators.get_diagonal_mask(an"CA")
    @test BitArray([false true true;true false true;true true false]) == gdm(pose).content

    @test [0.0 0.72614 0.68519; 0.72614 0.0 0.66099; 0.68519 0.66099 0.0] == ProtoSyn.Calculators.load_map("cmap-test-1.txt")

    # TODO
end