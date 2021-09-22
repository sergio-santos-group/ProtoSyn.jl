println("-----------\n Calculators:")

@testset verbose = true "Calculators $(repeat("-", 45))" begin

    @testset verbose = true "$(@sprintf "%-54s" "Verlet List")" begin
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

        vl4 = ProtoSyn.Calculators.VerletList(pose, an"CA"); vl4.cutoff = 8.0
        ProtoSyn.Calculators.update!(ProtoSyn.SISD_0, vl4, pose, an"CA")
        @test vl4.offset[1] === 1
        @test vl4.offset[3] === 6
        @test vl4.list[vl4.offset[3]] === -1

        vl5 = ProtoSyn.Calculators.VerletList(pose, an"CA"); vl5.cutoff = 8.0
        ProtoSyn.Calculators.update!(ProtoSyn.SIMD_1, vl5, pose, an"CA")

        @test vl4.offset == vl5.offset
        @test vl4.list[vl4.offset[3]] == vl5.list[vl5.offset[3]]
    end


    @testset verbose = true "$(@sprintf "%-54s" "SISD_0 Distance Matrix")" begin
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

        vl = ProtoSyn.Calculators.VerletList(pose); vl.cutoff = 8.0
        ProtoSyn.Calculators.update!(ProtoSyn.SISD_0, vl, pose)
        dm_vl1 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SISD_0, pose.state.x.coords, vl)
        dm_vl2 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SISD_0, pose.state, vl)
        dm_vl3 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SISD_0, pose, vl)
        @test dm_vl1 == dm_vl2
        @test dm_vl1 == dm_vl3
        @test dm_vl1[3, 10] ≈ 3.7194380493597423 atol = 1e-5
    end

    @testset verbose = true "$(@sprintf "%-54s" "SIMD_1 Distance Matrix")" begin
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

        vl = ProtoSyn.Calculators.VerletList(pose); vl.cutoff = 8.0
        ProtoSyn.Calculators.update!(ProtoSyn.SIMD_1, vl, pose)
        dm_vl1 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state.x.coords, vl)
        dm_vl2 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state, vl)
        dm_vl3 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose, vl)
        @test dm_vl1 == dm_vl2
        @test dm_vl1 == dm_vl3
        @test dm_vl1[3, 10] ≈ 3.7194380493597423 atol = 1e-5
    end

    if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2
        @testset verbose = true "$(@sprintf "%-54s" "CUDA_2 Distance Matrix")" begin
            pose = copy(backup)
            dm = collect(ProtoSyn.Calculators.distance_matrix(ProtoSyn.CUDA_2, pose.state.x.coords))

            @test dm[1, 1] === 0.0
            @test size(dm) === (39, 39)
            @test dm[end, 1] === dm[1, end]
            @test dm[1, end] ≈ 9.503650710901654 atol = 1e-5
            @test collect(ProtoSyn.Calculators.distance_matrix(ProtoSyn.CUDA_2, pose.state)) == dm
            @test collect(ProtoSyn.Calculators.distance_matrix(ProtoSyn.CUDA_2, pose)) == dm
            
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
            @test collect(ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.CUDA_2, pose.state)) == dm_full
            @test collect(ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.CUDA_2, pose)) == dm_full

            dm_full_sele = collect(ProtoSyn.Calculators.full_distance_matrix(ProtoSyn.CUDA_2, pose, an"CA"))
            @test dm_full_sele[1, 1] === 0.0
            @test size(dm_full_sele) === (3, 3)
            @test dm_full_sele[end, 1] === dm_full_sele[1, end]
            @test dm_full_sele[1, end] ≈ 7.065359511597963 atol = 1e-5

            vl = ProtoSyn.Calculators.VerletList(pose); vl.cutoff = 8.0
            ProtoSyn.Calculators.update!(ProtoSyn.SIMD_1, vl, pose)
            dm_vl1 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state.x.coords, vl)
            dm_vl2 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state, vl)
            dm_vl3 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose, vl)
            @test dm_vl1 == dm_vl2
            @test dm_vl1 == dm_vl3
            @test dm_vl1[3, 10] ≈ 3.7194380493597423 atol = 1e-5
        end
    else
        @warn "Skipping CUDA_2 Distance Matrix test on this system: No CUDA found."
    end

    @testset verbose = true "$(@sprintf "%-54s" "Dynamic-Acceleration-Type Distance Matrix")" begin
        pose = copy(backup)
        dm = collect(ProtoSyn.Calculators.distance_matrix(pose.state.x.coords))

        @test dm[1, 1] === 0.0
        @test size(dm) === (39, 39)
        @test dm[1, end] ≈ 9.503650710901654 atol = 1e-5
        @test collect(ProtoSyn.Calculators.distance_matrix(pose.state)) == dm
        @test collect(ProtoSyn.Calculators.distance_matrix(pose)) == dm
        
        dm_sele = collect(ProtoSyn.Calculators.distance_matrix(pose, an"CA"))
        @test dm_sele[1, 1] === 0.0
        @test size(dm_sele) === (3, 3)
        @test dm_sele[1, end] ≈ 7.065359511597963 atol = 1e-5

        dm_full = collect(ProtoSyn.Calculators.full_distance_matrix(pose.state.x.coords))
        @test dm_full[1, 1] === 0.0
        @test size(dm_full) === (39, 39)
        @test dm_full[1, end] ≈ 9.503650710901654 atol = 1e-5
        @test collect(ProtoSyn.Calculators.full_distance_matrix(pose.state)) == dm_full
        @test collect(ProtoSyn.Calculators.full_distance_matrix(pose)) == dm_full

        dm_full_sele = collect(ProtoSyn.Calculators.full_distance_matrix(pose, an"CA"))
        @test dm_full_sele[1, 1] === 0.0
        @test size(dm_full_sele) === (3, 3)
        @test dm_full_sele[1, end] ≈ 7.065359511597963 atol = 1e-5
    end

    @testset verbose = true "$(@sprintf "%-54s" "Potential")" begin
        
        # * DYNAMIC
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

        @test [0.0 0.72614 0.68519; 0.72614 0.0 0.66099; 0.68519 0.66099 0.0] == ProtoSyn.Calculators.load_map("Core/cmap-test-1.txt")

        p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.CUDA_2; d3=10.0, d4 = 20.0)
        @test p(100.0) === 1700.0
        e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, true, nothing)
        @test sum(e) ≈ 10.376791600496576 atol = 1e-10
        @test sum(f[:, 4, :]) ≈ -0.17472360105998336 atol = 1e-10

        e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, true, nothing)
        @test sum(e) ≈ 10.376791600496576 atol = 1e-10
        @test sum(f[:, 4, :]) ≈ -0.17472360105998336 atol = 1e-10

        p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.SIMD_1, d3=10.0, d4 = 20.0)
        @test p(100.0) === 1700.0
        e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, true, nothing)
        @test sum(e) ≈ 10.376791600496576 atol = 1e-10
        @test sum(f[:, 4, :]) ≈ -0.17472360105998336 atol = 1e-10

        # Mask function
        mask = ProtoSyn.Calculators.get_diagonal_mask()

        p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.CUDA_2, d3=10.0, d4 = 20.0)
        e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, true, nothing, mask(pose))
        @test e ≈ 10.376791600496576 atol = 1e-10
        @test sum(f[:, 4, :]) ≈ -0.17472360105998336 atol = 1e-10

        using Random
        Random.seed!(1)

        map = rand(pose.state.size, pose.state.size)

        @testset verbose = true "$(@sprintf "%-52s" "No Selection / Mask")" begin
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, true, nothing, map)
            @test e ≈ 5.086967385156593 atol = 1e-10
            @test f[1, 1] ≈ -1.8414104049560087 atol = 1e-10
            @test f[2, 1] ≈ 0.5557481135666688 atol = 1e-10
            @test f[3, 1] ≈ -0.5090176478479795 atol = 1e-10

            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, true, nothing, map)
            @test e ≈ 5.086967385156593 atol = 1e-10
            @test f[1, 1] ≈ -1.8414104049560087 atol = 1e-10
            @test f[2, 1] ≈ 0.5557481135666688 atol = 1e-10
            @test f[3, 1] ≈ -0.5090176478479795 atol = 1e-10

            p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.SIMD_1, d3=10.0, d4 = 20.0)
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, true, nothing, map)
            @test e ≈ 5.086967385156593 atol = 1e-10
            @test f[1, 1] ≈ -1.8414104049560087 atol = 1e-10
            @test f[2, 1] ≈ 0.5557481135666688 atol = 1e-10
            @test f[3, 1] ≈ -0.5090176478479795 atol = 1e-10

            function get_map(pose::Pose)
                Random.seed!(1)
                return rand(pose.state.size, pose.state.size)
            end

            @testset verbose = true "$(@sprintf "%-50s" "Mask == Function")" begin
                p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.CUDA_2, d3=10.0, d4 = 20.0)
                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, true, nothing, get_map)
                @test e ≈ 5.086967385156593 atol = 1e-10
                @test f[1, 1] ≈ -1.8414104049560087 atol = 1e-10
                @test f[2, 1] ≈ 0.5557481135666688 atol = 1e-10
                @test f[3, 1] ≈ -0.5090176478479795 atol = 1e-10

                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, true, nothing, get_map)
                @test e ≈ 5.086967385156593 atol = 1e-10
                @test f[1, 1] ≈ -1.8414104049560087 atol = 1e-10
                @test f[2, 1] ≈ 0.5557481135666688 atol = 1e-10
                @test f[3, 1] ≈ -0.5090176478479795 atol = 1e-10

                p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.SIMD_1, d3=10.0, d4 = 20.0)
                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, true, nothing, get_map)
                @test e ≈ 5.086967385156593 atol = 1e-10
                @test f[1, 1] ≈ -1.8414104049560087 atol = 1e-10
                @test f[2, 1] ≈ 0.5557481135666688 atol = 1e-10
                @test f[3, 1] ≈ -0.5090176478479795 atol = 1e-10
            end
        end

        @testset verbose = true "$(@sprintf "%-52s" "No Selection / No Mask")" begin
            p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.CUDA_2, d3=10.0, d4 = 20.0)
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, true, nothing)
            @test e ≈ 10.376791600496576 atol = 1e-10
            @test sum(f[:, 4, :]) ≈ -0.17472360105998336 atol = 1e-10

            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, true, nothing)
            @test e ≈ 10.376791600496576 atol = 1e-10
            @test sum(f[:, 4, :]) ≈ -0.17472360105998336 atol = 1e-10

            p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.SIMD_1, d3=10.0, d4 = 20.0)
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, true, nothing)
            @test e ≈ 10.376791600496576 atol = 1e-10
            @test sum(f[:, 4, :]) ≈ -0.17472360105998336 atol = 1e-10
        end

        @testset verbose = true "$(@sprintf "%-52s" "Selection / No Mask")" begin
            p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.CUDA_2, d1=2.0, d2 = 4.0)
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, true, nothing, an"CA")
            @test e ≈ 0.1573131450006226 atol = 1e-10
            @test f[1, 3] ≈ -0.5364771422229443 atol = 1e-10
            @test f[2, 3] ≈ 0.16304853046417908 atol = 1e-10
            @test f[3, 3] ≈ -1.973899021248059e-8 atol = 1e-10

            
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, true, nothing, an"CA")
            @test e ≈ 0.1573131450006226 atol = 1e-10
            @test f[1, 3] ≈ -0.5364771422229443 atol = 1e-10
            @test f[2, 3] ≈ 0.16304853046417908 atol = 1e-10
            @test f[3, 3] ≈ -1.973899021248059e-8 atol = 1e-10

            p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.SIMD_1, d1=2.0, d2 = 4.0)
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, true, nothing, an"CA")
            @test e ≈ 0.1573131450006226 atol = 1e-10
            @test f[1, 3] ≈ -0.5364771422229443 atol = 1e-10
            @test f[2, 3] ≈ 0.16304853046417908 atol = 1e-10
            @test f[3, 3] ≈ -1.973899021248059e-8 atol = 1e-10
        end

        mask = ProtoSyn.Calculators.diagonal_mask(pose, an"CA")

        @testset verbose = true "$(@sprintf "%-52s" "Selection / Mask")" begin
            p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.CUDA_2, d1=2.0, d2 = 4.0)
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, true, nothing, an"CA", mask)
            @test e ≈ 0.1573131450006226 atol = 1e-10
            @test f[1, 3] ≈ -0.5364771422229443 atol = 1e-10
            @test f[2, 3] ≈ 0.16304853046417908 atol = 1e-10
            @test f[3, 3] ≈ -1.973899021248059e-8 atol = 1e-10

            
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, true, nothing, an"CA", mask)
            @test e ≈ 0.1573131450006226 atol = 1e-10
            @test f[1, 3] ≈ -0.5364771422229443 atol = 1e-10
            @test f[2, 3] ≈ 0.16304853046417908 atol = 1e-10
            @test f[3, 3] ≈ -1.973899021248059e-8 atol = 1e-10

            p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.SIMD_1, d1=2.0, d2 = 4.0)
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, true, nothing, an"CA", mask)
            @test e ≈ 0.1573131450006226 atol = 1e-10
            @test f[1, 3] ≈ -0.5364771422229443 atol = 1e-10
            @test f[2, 3] ≈ 0.16304853046417908 atol = 1e-10
            @test f[3, 3] ≈ -1.973899021248059e-8 atol = 1e-10
        
            @testset verbose = true "$(@sprintf "%-50s" "Mask == Function")" begin

                mask = ProtoSyn.Calculators.get_diagonal_mask(an"CA")

                p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.CUDA_2, d1=2.0, d2 = 4.0)
                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, true, nothing, an"CA", mask)
                @test e ≈ 0.1573131450006226 atol = 1e-10
                @test f[1, 3] ≈ -0.5364771422229443 atol = 1e-10
                @test f[2, 3] ≈ 0.16304853046417908 atol = 1e-10
                @test f[3, 3] ≈ -1.973899021248059e-8 atol = 1e-10

                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, true, nothing, an"CA", mask)
                @test e ≈ 0.1573131450006226 atol = 1e-10
                @test f[1, 3] ≈ -0.5364771422229443 atol = 1e-10
                @test f[2, 3] ≈ 0.16304853046417908 atol = 1e-10
                @test f[3, 3] ≈ -1.973899021248059e-8 atol = 1e-10

                p = ProtoSyn.Calculators.get_flat_bottom_potential(ProtoSyn.SIMD_1, d1=2.0, d2 = 4.0)
                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, true, nothing, an"CA", mask)
                @test e ≈ 0.1573131450006226 atol = 1e-10
                @test f[1, 3] ≈ -0.5364771422229443 atol = 1e-10
                @test f[2, 3] ≈ 0.16304853046417908 atol = 1e-10
                @test f[3, 3] ≈ -1.973899021248059e-8 atol = 1e-10
            end
        end
    end

    @testset verbose = true "$(@sprintf "%-54s" "TorchANI")" begin
        pose = copy(backup)
        
        # get_ani_species
        @test ProtoSyn.Calculators.TorchANI.get_ani_species(pose.graph[1][1]) == [7, 1, 6, 1, 1, 6, 8]
        @test ProtoSyn.Calculators.TorchANI.get_ani_species(pose)[1:7] == [7, 1, 6, 1, 1, 6, 8]

        # --- model

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model(ProtoSyn.CUDA_2, pose)
        @test e ≈ -0.1257355809211731 atol = 1e-5
        @test f === nothing

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model(ProtoSyn.SIMD_1, pose)
        @test e ≈ -0.1257355809211731 atol = 1e-5
        @test f === nothing

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model(ProtoSyn.SISD_0, pose)
        @test e ≈ -0.1257355809211731 atol = 1e-5
        @test f === nothing

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model(ProtoSyn.CUDA_2, pose, true)
        @test f[1, 1] ≈ 0.10562094300985336 atol = 1e-5
        @test f[2, 1] ≈ 0.08480343967676163 atol = 1e-5
        @test f[3, 1] ≈ -0.0004017162136733532 atol = 1e-5

        m = ProtoSyn.Calculators.TorchANI.get_default_torchani_model()
        @test m.calc === ProtoSyn.Calculators.TorchANI.calc_torchani_model

        # --- ensemble

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(ProtoSyn.CUDA_2, pose)
        @test e ≈ -0.12801788747310638 atol = 1e-5
        @test f === nothing

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(ProtoSyn.SIMD_1, pose)
        @test e ≈ -0.12801788747310638 atol = 1e-5
        @test f === nothing

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(ProtoSyn.SISD_0, pose)
        @test e ≈ -0.12801788747310638 atol = 1e-5
        @test f === nothing

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(ProtoSyn.CUDA_2, pose, true)
        @test f[1, 1] ≈ 0.1022588387131691 atol = 1e-5
        @test f[2, 1] ≈ 0.1050674319267273 atol = 1e-5
        @test f[3, 1] ≈ -0.000596923753619194 atol = 1e-5

        m = ProtoSyn.Calculators.TorchANI.get_default_torchani_ensemble()
        @test m.calc === ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble

        # --- model xml-rpc

        test_xml_rpc = false # ! Setting to false because can't pass Travis CI
        if test_xml_rpc
            @test ProtoSyn.Calculators.TorchANI.server === nothing
            ProtoSyn.Calculators.TorchANI.start_torchANI_server()
            @test typeof(ProtoSyn.Calculators.TorchANI.server) === Base.Process
            ProtoSyn.Calculators.TorchANI.stop_torchANI_server()
            @test ProtoSyn.Calculators.TorchANI.server === nothing

            @test_throws ErrorException ProtoSyn.Calculators.TorchANI.calc_torchani_model_xmlrpc(ProtoSyn.SISD_0, pose)
            @test_throws ErrorException ProtoSyn.Calculators.TorchANI.calc_torchani_model_xmlrpc(ProtoSyn.SIMD_1, pose)

            e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model_xmlrpc(ProtoSyn.CUDA_2, pose)
            @test e ≈ -0.12573562562465668 atol = 1e-5
            @test f === nothing

            if ProtoSyn.acceleration.active == ProtoSyn.CUDA_2
                e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model_xmlrpc(ProtoSyn.CUDA_2, pose)
                @test e ≈ -0.12573561072349548 atol = 1e-5
                @test f === nothing
            end

            e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model_xmlrpc(ProtoSyn.CUDA_2, pose, true)
            @test f[1, 1] ≈ 0.10562092810869217 atol = 1e-5
            @test f[2, 1] ≈ 0.08480347692966461 atol = 1e-5
            @test f[3, 1] ≈ -0.0004017228784505278 atol = 1e-5
            
            if ProtoSyn.acceleration.active == ProtoSyn.CUDA_2
                e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model_xmlrpc(ProtoSyn.CUDA_2, pose, true)
                @test f[1, 1] ≈ 0.10562092810869217 atol = 1e-5
                @test f[2, 1] ≈ 0.08480347692966461 atol = 1e-5
                @test f[3, 1] ≈ -0.0004017228784505278 atol = 1e-5
            end

            m = ProtoSyn.Calculators.TorchANI.get_default_torchani_model_xmlrpc()
            @test m.calc === ProtoSyn.Calculators.TorchANI.calc_torchani_model_xmlrpc
        end
    end
end