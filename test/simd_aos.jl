@generated function simd_aos(state::State{T}, vlist::VerletList, ::Type{Val{Forces}}, ::Type{Vec{4,T}}) where {Forces,T<:AbstractFloat}
@assert isa(Forces, Bool) "Forces must be a boolean"
quote

    σ = one(T)
    ϵ = ones(T, 4)
    energy = zero(T)
    natoms = state.size

    mlane = VecRange{4}(1)  # mask lane
    coords = view(state.coords, :)
    cutsq = convert(T, vlist.cutoff*vlist.cutoff)
    
    #region FORCE_SECTION
    $(if Forces === true
        quote
            flane = VecRange{4}(0)  # force lane
            forces = view(state.forces, :)
            # fi = zeros(T, 4)
        end
    end)
    #endregion FORCE_SECTION

    m1 = (0, 4, 1, 5)
    m2 = (2, 6, 3, 7)
    mx = (0, 1, 4, 5)
    my = (2, 3, 6, 7)

    @inbounds for i = 1:natoms-1
        
        ptr = vlist.offset[i]
        if vlist.list[ptr] < 1
            continue
        end
        ptr_stop = vlist.offset[i+1]-2

        $(if Forces=== true
            :(fi = Vec{4,T}(0))
            # :(fill!(fi, T(0)))
        end)

        i1 = i<<2 - 3
        vi = vload(Vec{4,T}, coords, i1)

        while ptr+3 <= ptr_stop
            j = ptr
            ptr += 4

            # j = 4*(j-1) + 1 = 4j-4+1 = 4j-3 = j<<2 - 3
            j_1 = vlist.list[j  ]<<2 - 3
            j_2 = vlist.list[j+1]<<2 - 3
            j_3 = vlist.list[j+2]<<2 - 3
            j_4 = vlist.list[j+3]<<2 - 3

            # println("$i  $(vlist.list[j])  $(vlist.list[j+1])  $(vlist.list[j+2])  $(vlist.list[j+3])")

            v1 = vload(Vec{4,T}, coords, j_1)
            v2 = vload(Vec{4,T}, coords, j_2)
            v3 = vload(Vec{4,T}, coords, j_3)
            v4 = vload(Vec{4,T}, coords, j_4)
            
            vi_1 = v1 - vi # xi1, yi1, zi1, wi1
            vi_2 = v2 - vi # xi2, yi2, zi2, wi2
            vi_3 = v3 - vi # xi3, yi3, zi3, wi3
            vi_4 = v4 - vi # xi4, yi4, zi4, wi4

            vs12_m1 = shufflevector(vi_1, vi_2, Val{m1})    # xi1, xi2, yi1, yi2
            vs34_m1 = shufflevector(vi_3, vi_4, Val{m1})    # xi3, xi4, yi3, yi4
            vs12_m2 = shufflevector(vi_1, vi_2, Val{m2})    # zi1, zi2, ?i1, ?i2
            vs34_m2 = shufflevector(vi_3, vi_4, Val{m2})    # zi3, zi4, ?i3, ?i4
            
            xij = shufflevector(vs12_m1, vs34_m1, Val{mx})   # xi1, xi2, xi3, xi4
            yij = shufflevector(vs12_m1, vs34_m1, Val{my})   # yi1, yi2, yi3, yi4
            zij = shufflevector(vs12_m2, vs34_m2, Val{mx})   # zi1, zi2, zi3, zi4

            dij_sq = xij*xij + yij*yij + zij*zij
            mask = dij_sq <= cutsq
            !any(mask) && continue

            ϵm = ϵ[mlane, mask]

            lj2 = σ / dij_sq
            lj6 = lj2*lj2*lj2
            energy += sum( ϵm*(lj6*lj6 - lj6) )
            
            #region FORCE_SECTION
            $(if Forces === true
                quote
                    fc = ϵm * T(24.0) * (lj6 - T(2.0) * lj6 * lj6) / dij_sq
                    @nexprs 4 u ->  if mask[u]
                        f = fc[u] * vi_u
                        forces[flane+j_u] -= f
                        fi += f
                    end
                    # if mask[1]
                    #     f = fc[1] * vi1
                    #     fi += f
                    #     forces[flane+j1] -= f
                    # end
                    # if mask[2]
                    #     f = fc[2] * vi2
                    #     fi += f
                    #     forces[flane+j2] -= f
                    # end
                    # if mask[3]
                    #     f = fc[3] * vi3
                    #     fi += f
                    #     forces[flane+j3] -= f
                    # end
                    # if mask[4]
                    #     f = fc[4] * vi4
                    #     fi += f
                    #     forces[flane+j4] -= f
                    # end
                end
            end)
            #endregion FORCE_SECTION

        end # while ptr
    
        # do remaining pairs
        for j = ptr:ptr_stop
            # println("$i  $(vlist.list[j])")
            j1 = vlist.list[j]<<2 - 3
            v1 = vload(Vec{4,T}, coords, j1)
            vi1 = v1 - vi
            
            sdij_sq = sum(vi1*vi1)
            (sdij_sq > cutsq) && continue

            slj2 = σ / sdij_sq
            slj6 = slj2*slj2*slj2
            energy += slj6*slj6 - slj6

            $(if Forces === true
                quote
                    sfc = T(24.0) * (slj6 - T(2.0) * slj6 * slj6) / sdij_sq
                    f = sfc * vi1
                    forces[flane + j1] -= f
                    fi += f
                end
            end)
        end
        # println("")

        #region FORCE_SECTION
        $(if Forces === true
            :(forces[flane + i1] += fi)
        end)
        #endregion FORCE_SECTION

    end # for i

end # quote
end # generated function














@generated function simd_aos(state::State{T}, vlist::VerletList, ::Type{Val{Forces}}, ::Type{Vec{8,T}}) where {Forces,T<:AbstractFloat}
    @assert isa(Forces, Bool) "Forces must be a boolean"
    quote
    
        σ = T(1)
        ϵ = ones(T, 8)
        energy = T(0)
        natoms = state.size
    
        mlane = VecRange{8}(1)  # mask lane
        coords = view(state.coords, :)
        cutsq = convert(T, vlist.cutoff*vlist.cutoff)
        
        #region FORCE_SECTION
        $(if Forces === true
            quote
                flane  = VecRange{4}(0)  # force lane
                forces = view(state.forces, :)
            end
        end)
        #endregion FORCE_SECTION
    
        m1 = (0, 4,  8, 12,  1,  5,  9, 13)
        m2 = (2, 6, 10, 14,  3,  7, 11, 15)
        mx = (0, 1,  2,  3,  8,  9, 10, 11)
        my = (4, 5,  6,  7, 12, 13, 14, 15)
        mj = (0,1,2,3,4,5,6,7)
        # mj = (0,1,2,3,8,9,10,11)
        
        @inbounds for i = 1:natoms-1
            
            ptr = vlist.offset[i]
            if vlist.list[ptr] < 1
                continue
            end
            ptr_stop = vlist.offset[i+1]-2
    
            $(if Forces=== true
                # :(fill!(fi, T(0)))
                # :(fi = Vec{8,T}(0))
                :(fi = Vec{4,T}(0))
            end)
    
            i1 = i<<2 - 3
            # vi = vload(Vec{8,T}, coords, i1)
            # vi = vgather(coords, Vec{8,Int}((0,1,2,3,0,1,2,3)) + i1)
            vi = vload(Vec{4,T}, coords, i1)

            while ptr+7 <= ptr_stop
                j = ptr
                ptr += 8
    
                # j = 4*(j-1) + 1 = 4j-4+1 = 4j-3 = j<<2 - 3
                
                j_1 = vlist.list[j  ]<<2 - 3
                j_2 = vlist.list[j+1]<<2 - 3
                j_3 = vlist.list[j+2]<<2 - 3
                j_4 = vlist.list[j+3]<<2 - 3
                j_5 = vlist.list[j+4]<<2 - 3
                j_6 = vlist.list[j+5]<<2 - 3
                j_7 = vlist.list[j+6]<<2 - 3
                j_8 = vlist.list[j+7]<<2 - 3

                vi_1 = vload(Vec{4,T}, coords, j_1) - vi
                vi_2 = vload(Vec{4,T}, coords, j_2) - vi
                vi_3 = vload(Vec{4,T}, coords, j_3) - vi
                vi_4 = vload(Vec{4,T}, coords, j_4) - vi
                vi_5 = vload(Vec{4,T}, coords, j_5) - vi
                vi_6 = vload(Vec{4,T}, coords, j_6) - vi
                vi_7 = vload(Vec{4,T}, coords, j_7) - vi
                vi_8 = vload(Vec{4,T}, coords, j_8) - vi
                
                vi12 = shufflevector(vi_1, vi_2, Val{mj})
                vi34 = shufflevector(vi_3, vi_4, Val{mj})
                vi56 = shufflevector(vi_5, vi_6, Val{mj})
                vi78 = shufflevector(vi_7, vi_8, Val{mj})

                #println("$i  $(vlist.list[j]) : $(vlist.list[j+7])")
    
                # v12 = vload(Vec{8,T}, coords, j12)
                # v34 = vload(Vec{8,T}, coords, j34)
                # v56 = vload(Vec{8,T}, coords, j56)
                # v78 = vload(Vec{8,T}, coords, j78)
                
                #vi12 = v12 - vi # xi1, yi1, zi1, wi1, xi2, yi2, zi2, wi2
                #vi34 = v34 - vi # xi3, yi3, zi3, wi3, xi4, yi4, zi4, wi4
                #vi56 = v56 - vi # xi5, yi5, zi5, wi5, xi6, yi6, zi6, wi6
                #vi78 = v78 - vi # xi7, yi7, zi7, wi7, xi8, yi8, zi8, wi8
                
                vs1234_m1 = shufflevector(vi12, vi34, Val{m1})    # xi1, xi2, xi3, xi4, yi1, yi2, yi3, yi4
                vs5678_m1 = shufflevector(vi56, vi78, Val{m1})    # xi5, xi6, xi7, xi8, yi5, yi6, yi7, yi8
                
                vs1234_m2 = shufflevector(vi12, vi34, Val{m1})    # zi1, zi2, zi3, zi4, wi1, wi2, wi3, wi4
                vs5678_m2 = shufflevector(vi56, vi78, Val{m1})    # zi5, zi6, zi7, zi8, wi5, wi6, wi7, wi8

                xij = shufflevector(vs1234_m1, vs5678_m1, Val{mx})   # xi1, xi2, xi3, xi4, xi5, xi6, xi7, xi8
                yij = shufflevector(vs1234_m1, vs5678_m1, Val{my})   # yi1, yi2, yi3, yi4, yi5, yi6, yi7, yi8
                zij = shufflevector(vs1234_m2, vs5678_m2, Val{mx})   # zi1, zi2, zi3, zi4, zi5, zi6, zi7, zi8

                dij_sq = xij*xij + yij*yij + zij*zij
                mask = dij_sq <= cutsq
                !any(mask) && continue
    
                
                lj2 = σ / dij_sq
                lj6 = lj2*lj2*lj2
                ϵm = ϵ[mlane, mask]
                energy += sum( ϵm*(lj6*lj6 - lj6) )
                
                # region FORCE_SECTION
                $(if Forces === true
                    quote
                        fc = ϵm * T(24.0) * (lj6 - T(2.0) * lj6 * lj6) / dij_sq
                        @nexprs 4 u ->  if mask[u]
                            f = fc[u] * vi_u
                            forces[flane+j_u] -= f
                            fi += f
                        end

                        # # fc = T(24.0) * (lj6 - T(2.0) * lj6 * lj6) / dij_sq
                        # #if mask[1]
                        #     # f1 = fc[1] * vi1
                        #     fi += fc[1] * vi1
                        #     forces[flane+j1] -= fc[1] * vi1
                        # #end
                        # #if mask[2]
                        #     # f2 = fc[2] * vi2
                        #     fi += fc[2] * vi2
                        #     forces[flane+j2] -= fc[2] * vi2
                        # #end
                        # #if mask[3]
                        #     # f3 = fc[3] * vi3
                        #     fi += fc[3] * vi3
                        #     forces[flane+j3] -= fc[3] * vi3
                        # #end
                        # #if mask[4]
                        #     # f4 = fc[4] * vi4
                        #     fi += fc[4] * vi4
                        #     forces[flane+j4] -= fc[4] * vi4
                        # #end
                        # #if mask[5]
                        #     # f5 = fc[5] * vi5
                        #     fi += fc[5] * vi5
                        #     forces[flane+j5] -= fc[5] * vi5
                        # #end
                        # #if mask[6]
                        #     # f6 = fc[6] * vi6
                        #     fi += fc[6] * vi6
                        #     forces[flane+j6] -= fc[6] * vi6
                        # #end
                        # #if mask[7]
                        #     # f7 = fc[7] * vi7
                        #     fi += fc[7] * vi7
                        #     forces[flane+j7] -= fc[7] * vi7
                        # #end
                        # #if mask[8]
                        #     # f8 = fc[8] * vi8
                        #     fi += fc[8] * vi8
                        #     forces[flane+j8] -= fc[8] * vi8
                        # #end
                    end
                end)
                #endregion FORCE_SECTION
    
            end # for i
            

            # do remaining pairs (atom-by-atom)
            vii = vload(Vec{4,T}, coords, i1)
            for j = ptr:ptr_stop
                #println("$i  $(vlist.list[j])")
                j1 = vlist.list[j]<<2 - 3
                vi1 = vload(Vec{4,T}, coords, j1) - vi

                sdij_sq = sum(vi1*vi1)
                (sdij_sq > cutsq) && continue
            
                slj2 = σ / sdij_sq
                slj6 = slj2*slj2*slj2
                energy += slj6*slj6 - slj6
                # energy += sum( slj6*slj6 - slj6 )
            
                $(if Forces === true
                    quote
                        sfc = T(24.0) * (slj6 - T(2.0) * slj6 * slj6) / sdij_sq
                        f = sfc * vi1
                        forces[flane + j1] -= f
                        fi += f
                    end
                end)
            end
            #println("")

            #region FORCE_SECTION
            $(if Forces === true
                :(forces[flane + i1] += fi)
            end)
            #endregion FORCE_SECTION
            
        end # for i
    
    end # quote
    end # generated function