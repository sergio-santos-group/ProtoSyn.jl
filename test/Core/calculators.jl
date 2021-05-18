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

        vl = ProtoSyn.Calculators.VerletList(pose.state.size); vl.cutoff = 8.0
        ProtoSyn.Calculators.update!(ProtoSyn.SISD_0, vl, pose, an"CA")
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

        vl = ProtoSyn.Calculators.VerletList(pose.state.size); vl.cutoff = 8.0
        ProtoSyn.Calculators.update!(ProtoSyn.SIMD_1, vl, pose, an"CA")
        dm_vl1 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state.x.coords, vl)
        dm_vl2 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose.state, vl)
        dm_vl3 = ProtoSyn.Calculators.distance_matrix(ProtoSyn.SIMD_1, pose, vl)
        @test dm_vl1 == dm_vl2
        @test dm_vl1 == dm_vl3
        @test dm_vl1[3, 10] ≈ 3.7194380493597423 atol = 1e-5
    end

    @testset verbose = true "$(@sprintf "%-54s" "CUDA_2 Distance Matrix")" begin
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

    @testset verbose = true "$(@sprintf "%-54s" "Dynamic-Acceleration-Type Distance Matrix")" begin
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

    @testset verbose = true "$(@sprintf "%-54s" "Potential")" begin
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

        e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose.state.x.coords[:], p)
        @test sum(e) ≈ 20.753583200992992 atol = 1e-10
        @test sum(f[:, 4, :]) ≈ 1.7650274195512723 atol = 1e-10

        @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose.state.x.coords[:], p)
        @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose.state.x.coords[:], p)

        if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
            e, f = ProtoSyn.Calculators.apply_potential(pose.state.x.coords[:], p)
            @test sum(e) ≈ 20.753583200992992 atol = 1e-10
            @test sum(f[:, 4, :]) ≈ 1.7650274195512723 atol = 1e-10
        end

        mask = ProtoSyn.Calculators.get_diagonal_mask(an"CA")

        @test_throws AssertionError ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose.state.x.coords[:], p, mask(pose))

        subset_pose = pose.state.x[:, findall(an"CA"(pose).content)]
        p = ProtoSyn.Calculators.get_flat_bottom_potential(d3=4.0, d4 = 5.0)
        e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, subset_pose[:], p, mask(pose))
        @test e ≈ 10.261438046391852 atol = 1e-10
        @test f[1, 1] ≈ 11.554902838944216 atol = 1e-10
        @test f[2, 1] ≈ -8.133968311658052 atol = 1e-10
        @test f[3, 1] ≈ -3.5777523521939767e-6 atol = 1e-10

        @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, subset_pose[:], p, mask(pose))
        @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, subset_pose[:], p, mask(pose))

        if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
            e, f = ProtoSyn.Calculators.apply_potential(subset_pose[:], p, mask(pose))
            @test e ≈ 10.261438046391852 atol = 1e-10
            @test f[1, 1] ≈ 11.554902838944216 atol = 1e-10
            @test f[2, 1] ≈ -8.133968311658052 atol = 1e-10
            @test f[3, 1] ≈ -3.5777523521939767e-6 atol = 1e-10
        end

        using Random
        Random.seed!(1)

        map = rand(pose.state.size, pose.state.size)

        @testset verbose = true "$(@sprintf "%-52s" "No Selection / Mask")" begin
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, map)
            @test e ≈ 1971.5700803399764 atol = 1e-10
            @test f[1, 1] ≈ 149.78398654397617 atol = 1e-10
            @test f[2, 1] ≈ -105.79414884147755 atol = 1e-10
            @test f[3, 1] ≈ -3.16249805532396 atol = 1e-10

            @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, map)
            @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, map)

            if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
                e, f = ProtoSyn.Calculators.apply_potential(pose, p, map)
                @test e ≈ 1971.5700803399764 atol = 1e-10
                @test f[1, 1] ≈ 149.78398654397617 atol = 1e-10
                @test f[2, 1] ≈ -105.79414884147755 atol = 1e-10
                @test f[3, 1] ≈ -3.16249805532396 atol = 1e-10
            end
        
            @testset verbose = true "$(@sprintf "%-50s" "Explicit selection == nothing")" begin
                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, nothing, map)
                @test e ≈ 1971.5700803399764 atol = 1e-10
                @test f[1, 1] ≈ 149.78398654397617 atol = 1e-10
                @test f[2, 1] ≈ -105.79414884147755 atol = 1e-10
                @test f[3, 1] ≈ -3.16249805532396 atol = 1e-10

                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, nothing, map)
                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, nothing, map)

                if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
                    e, f = ProtoSyn.Calculators.apply_potential(pose, p, nothing, map)
                    @test e ≈ 1971.5700803399764 atol = 1e-10
                    @test f[1, 1] ≈ 149.78398654397617 atol = 1e-10
                    @test f[2, 1] ≈ -105.79414884147755 atol = 1e-10
                    @test f[3, 1] ≈ -3.16249805532396 atol = 1e-10
                end
            end

            function get_map(pose::Pose)
                Random.seed!(1)
                return rand(pose.state.size, pose.state.size)
            end

            @testset verbose = true "$(@sprintf "%-50s" "Mask == Function")" begin
                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, get_map)
                @test e ≈ 1971.5700803399764 atol = 1e-10
                @test f[1, 1] ≈ 149.78398654397617 atol = 1e-10
                @test f[2, 1] ≈ -105.79414884147755 atol = 1e-10
                @test f[3, 1] ≈ -3.16249805532396 atol = 1e-10

                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, get_map)
                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, get_map)

                if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
                    e, f = ProtoSyn.Calculators.apply_potential(pose, p, get_map)
                    @test e ≈ 1971.5700803399764 atol = 1e-10
                    @test f[1, 1] ≈ 149.78398654397617 atol = 1e-10
                    @test f[2, 1] ≈ -105.79414884147755 atol = 1e-10
                    @test f[3, 1] ≈ -3.16249805532396 atol = 1e-10
                end
            end

            @testset verbose = true "$(@sprintf "%-50s" "Explicit selection == nothing && Mask == Function")" begin
                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, nothing, get_map)
                @test e ≈ 1971.5700803399764 atol = 1e-10
                @test f[1, 1] ≈ 149.78398654397617 atol = 1e-10
                @test f[2, 1] ≈ -105.79414884147755 atol = 1e-10
                @test f[3, 1] ≈ -3.16249805532396 atol = 1e-10

                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, nothing, get_map)
                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, nothing, get_map)

                if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
                    e, f = ProtoSyn.Calculators.apply_potential(pose, p, nothing, get_map)
                    @test e ≈ 1971.5700803399764 atol = 1e-10
                    @test f[1, 1] ≈ 149.78398654397617 atol = 1e-10
                    @test f[2, 1] ≈ -105.79414884147755 atol = 1e-10
                    @test f[3, 1] ≈ -3.16249805532396 atol = 1e-10
                end
            end
        end

        @testset verbose = true "$(@sprintf "%-52s" "No Selection / No Mask")" begin
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p)
            @test e ≈ 3961.850868099287 atol = 1e-10
            @test f[1, 1] ≈ 355.5495414738819 atol = 1e-10
            @test f[2, 1] ≈ -258.4799542670788 atol = 1e-10
            @test f[3, 1] ≈ -6.149467896237443 atol = 1e-10

            @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p)
            @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p)

            if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
                e, f = ProtoSyn.Calculators.apply_potential(pose, p)
                @test e ≈ 3961.850868099287 atol = 1e-10
                @test f[1, 1] ≈ 355.5495414738819 atol = 1e-10
                @test f[2, 1] ≈ -258.4799542670788 atol = 1e-10
                @test f[3, 1] ≈ -6.149467896237443 atol = 1e-10
            end
        
            @testset verbose = true "$(@sprintf "%-50s" "Explicit selection == nothing")" begin
                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, nothing)
                @test e ≈ 3961.850868099287 atol = 1e-10
                @test f[1, 1] ≈ 355.5495414738819 atol = 1e-10
                @test f[2, 1] ≈ -258.4799542670788 atol = 1e-10
                @test f[3, 1] ≈ -6.149467896237443 atol = 1e-10

                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, nothing)
                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, nothing)

                if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
                    e, f = ProtoSyn.Calculators.apply_potential(pose, p, nothing)
                    @test e ≈ 3961.850868099287 atol = 1e-10
                    @test f[1, 1] ≈ 355.5495414738819 atol = 1e-10
                    @test f[2, 1] ≈ -258.4799542670788 atol = 1e-10
                    @test f[3, 1] ≈ -6.149467896237443 atol = 1e-10
                end
            end

            @testset verbose = true "$(@sprintf "%-50s" "Explicit selection == nothing && mask == nothing")" begin
                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, nothing, nothing)
                @test e ≈ 3961.850868099287 atol = 1e-10
                @test f[1, 1] ≈ 355.5495414738819 atol = 1e-10
                @test f[2, 1] ≈ -258.4799542670788 atol = 1e-10
                @test f[3, 1] ≈ -6.149467896237443 atol = 1e-10

                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, nothing, nothing)
                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, nothing, nothing)

                if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
                    e, f = ProtoSyn.Calculators.apply_potential(pose, p, nothing, nothing)
                    @test e ≈ 3961.850868099287 atol = 1e-10
                    @test f[1, 1] ≈ 355.5495414738819 atol = 1e-10
                    @test f[2, 1] ≈ -258.4799542670788 atol = 1e-10
                    @test f[3, 1] ≈ -6.149467896237443 atol = 1e-10
                end
            end
        end

        @testset verbose = true "$(@sprintf "%-52s" "Selection / No Mask")" begin
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, an"CA")
            @test e ≈ 10.261438046391852 atol = 1e-10
            @test f[1, 3] ≈ 11.554902838944216 atol = 1e-10
            @test f[2, 3] ≈ -8.133968311658052 atol = 1e-10
            @test f[3, 3] ≈ -3.5777523521939767e-6 atol = 1e-10

            @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, an"CA")
            @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, an"CA")

            if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
                e, f = ProtoSyn.Calculators.apply_potential(pose, p, an"CA")
                @test e ≈ 10.261438046391852 atol = 1e-10
                @test f[1, 3] ≈ 11.554902838944216 atol = 1e-10
                @test f[2, 3] ≈ -8.133968311658052 atol = 1e-10
                @test f[3, 3] ≈ -3.5777523521939767e-6 atol = 1e-10
            end
        
            @testset verbose = true "$(@sprintf "%-50s" "Explicit mask == nothing")" begin
                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, an"CA", nothing)
                @test e ≈ 10.261438046391852 atol = 1e-10
                @test f[1, 3] ≈ 11.554902838944216 atol = 1e-10
                @test f[2, 3] ≈ -8.133968311658052 atol = 1e-10
                @test f[3, 3] ≈ -3.5777523521939767e-6 atol = 1e-10

                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, an"CA", nothing)
                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, an"CA", nothing)

                if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
                    e, f = ProtoSyn.Calculators.apply_potential(pose, p, an"CA", nothing)
                    @test e ≈ 10.261438046391852 atol = 1e-10
                    @test f[1, 3] ≈ 11.554902838944216 atol = 1e-10
                    @test f[2, 3] ≈ -8.133968311658052 atol = 1e-10
                    @test f[3, 3] ≈ -3.5777523521939767e-6 atol = 1e-10
                end
            end
        end

        mask = ProtoSyn.Calculators.diagonal_mask(pose, an"CA")

        @testset verbose = true "$(@sprintf "%-52s" "Selection / Mask")" begin
            e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, an"CA", mask)
            @test e ≈ 10.261438046391852 atol = 1e-10
            @test f[1, 3] ≈ 11.554902838944216 atol = 1e-10
            @test f[2, 3] ≈ -8.133968311658052 atol = 1e-10
            @test f[3, 3] ≈ -3.5777523521939767e-6 atol = 1e-10

            @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, an"CA", mask)
            @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, an"CA", mask)

            if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
                e, f = ProtoSyn.Calculators.apply_potential(pose, p, an"CA", mask)
                @test e ≈ 10.261438046391852 atol = 1e-10
                @test f[1, 3] ≈ 11.554902838944216 atol = 1e-10
                @test f[2, 3] ≈ -8.133968311658052 atol = 1e-10
                @test f[3, 3] ≈ -3.5777523521939767e-6 atol = 1e-10
            end
        
            @testset verbose = true "$(@sprintf "%-50s" "Mask == Function")" begin

                mask = ProtoSyn.Calculators.get_diagonal_mask(an"CA")

                e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, p, an"CA", mask)
                @test e ≈ 10.261438046391852 atol = 1e-10
                @test f[1, 3] ≈ 11.554902838944216 atol = 1e-10
                @test f[2, 3] ≈ -8.133968311658052 atol = 1e-10
                @test f[3, 3] ≈ -3.5777523521939767e-6 atol = 1e-10

                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SISD_0, pose, p, an"CA", mask)
                @test_throws ErrorException ProtoSyn.Calculators.apply_potential(ProtoSyn.SIMD_1, pose, p, an"CA", mask)

                if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2 # * DYNAMIC
                    e, f = ProtoSyn.Calculators.apply_potential(pose, p, an"CA", mask)
                    @test e ≈ 10.261438046391852 atol = 1e-10
                    @test f[1, 3] ≈ 11.554902838944216 atol = 1e-10
                    @test f[2, 3] ≈ -8.133968311658052 atol = 1e-10
                    @test f[3, 3] ≈ -3.5777523521939767e-6 atol = 1e-10
                end
            end
        end
    end

    @testset verbose = true "$(@sprintf "%-54s" "TorchANI")" begin
        pose = copy(backup)
        
        # get_ani_species
        @test ProtoSyn.Calculators.TorchANI.get_ani_species(pose.graph[1][1]) == [7, 1, 6, 1, 1, 6, 8]
        @test ProtoSyn.Calculators.TorchANI.get_ani_species(pose)[1:7] == [7, 1, 6, 1, 1, 6, 8]

        # --- model

        @test_throws ErrorException ProtoSyn.Calculators.TorchANI.calc_torchani_model(ProtoSyn.SISD_0, pose)
        @test_throws ErrorException ProtoSyn.Calculators.TorchANI.calc_torchani_model(ProtoSyn.SIMD_1, pose)

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model(ProtoSyn.CUDA_2, pose)
        @test e ≈ -0.1257355809211731 atol = 1e-5
        @test f === nothing

        if ProtoSyn.acceleration.active == ProtoSyn.CUDA_2 # * DYNAMIC
            e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model(pose)
            @test e ≈ -0.1257355809211731 atol = 1e-5
            @test f === nothing
        end

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model(ProtoSyn.CUDA_2, pose, true)
        @test f[1, 1] ≈ 0.10562094300985336 atol = 1e-5
        @test f[2, 1] ≈ 0.08480343967676163 atol = 1e-5
        @test f[3, 1] ≈ -0.0004017162136733532 atol = 1e-5

        if ProtoSyn.acceleration.active == ProtoSyn.CUDA_2 # * DYNAMIC
            e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_model(pose, true)
            @test f[1, 1] ≈ 0.10562094300985336 atol = 1e-5
            @test f[2, 1] ≈ 0.08480343967676163 atol = 1e-5
            @test f[3, 1] ≈ -0.0004017162136733532 atol = 1e-5
        end

        m = ProtoSyn.Calculators.TorchANI.get_default_torchani_model()
        @test m.calc === ProtoSyn.Calculators.TorchANI.calc_torchani_model

        # --- ensemble

        @test_throws ErrorException ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(ProtoSyn.SISD_0, pose)
        @test_throws ErrorException ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(ProtoSyn.SIMD_1, pose)

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(ProtoSyn.CUDA_2, pose)
        @test e ≈ -0.12801788747310638 atol = 1e-5
        @test f === nothing

        if ProtoSyn.acceleration.active == ProtoSyn.CUDA_2
            e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(ProtoSyn.CUDA_2, pose)
            @test e ≈ -0.12801790237426758 atol = 1e-5
            @test f === nothing
        end

        e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(ProtoSyn.CUDA_2, pose, true)
        @test f[1, 1] ≈ 0.1022588387131691 atol = 1e-5
        @test f[2, 1] ≈ 0.1050674319267273 atol = 1e-5
        @test f[3, 1] ≈ -0.000596923753619194 atol = 1e-5
        
        if ProtoSyn.acceleration.active == ProtoSyn.CUDA_2
            e, f = ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(ProtoSyn.CUDA_2, pose, true)
            @test f[1, 1] ≈ 0.1022588387131691 atol = 1e-5
            @test f[2, 1] ≈ 0.1050674319267273 atol = 1e-5
            @test f[3, 1] ≈ -0.000596923753619194 atol = 1e-5
        end

        m = ProtoSyn.Calculators.TorchANI.get_default_torchani_ensemble()
        @test m.calc === ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble

        # --- model xml-rpc

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