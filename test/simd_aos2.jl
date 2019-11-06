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
        end)

        i1 = i<<2 - 3
        vi = vload(Vec{4,T}, coords, i1)

        while ptr+3 <= ptr_stop
            j = ptr
            ptr += 4

            # j = 4*(j-1) + 1 = 4j-4+1 = 4j-3 = j<<2 - 3
            j1 = vlist.list[j  ]<<2 - 3
            j2 = vlist.list[j+1]<<2 - 3
            j3 = vlist.list[j+2]<<2 - 3
            j4 = vlist.list[j+3]<<2 - 3
            # println("$i  $(vlist.list[j])  $(vlist.list[j+1])  $(vlist.list[j+2])  $(vlist.list[j+3])")

            v1 = vload(Vec{4,T}, coords, j1)
            v2 = vload(Vec{4,T}, coords, j2)
            v3 = vload(Vec{4,T}, coords, j3)
            v4 = vload(Vec{4,T}, coords, j4)
            
            vi1 = v1 - vi # xi1, yi1, zi1, wi1
            vi2 = v2 - vi # xi2, yi2, zi2, wi2
            vi3 = v3 - vi # xi3, yi3, zi3, wi3
            vi4 = v4 - vi # xi4, yi4, zi4, wi4

            vs12_m1 = shufflevector(vi1, vi2, Val{m1})    # xi1, xi2, yi1, yi2
            vs34_m1 = shufflevector(vi3, vi4, Val{m1})    # xi3, xi4, yi3, yi4
            vs12_m2 = shufflevector(vi1, vi2, Val{m2})    # zi1, zi2, ?i1, ?i2
            vs34_m2 = shufflevector(vi3, vi4, Val{m2})    # zi3, zi4, ?i3, ?i4
            
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
                    #if mask[1]
                        fi += fc[1] * vi1
                        forces[flane+j1] -= fc[1] * vi1
                    #end
                    #if mask[2]
                        fi += fc[2] * vi2
                        forces[flane+j2] -= fc[2] * vi2
                    #end
                    #if mask[3]
                        fi += fc[3] * vi3
                        forces[flane+j3] -= fc[3] * vi3
                    #end
                    #if mask[4]
                        fi += fc[4] * vi4
                        forces[flane+j4] -= fc[4] * vi4
                    #end
                    # MUST CHECK CUTOFF MASK!!!!!!
                end
            end)
            #endregion FORCE_SECTION

        end # for i
    
        #region FORCE_SECTION
        $(if Forces === true
            :(forces[flane + i1] += fi)
        end)
        #endregion FORCE_SECTION

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
            # energy += sum( slj6*slj6 - slj6 )
            # println(typeof(slj2))

            $(if Forces === true
                quote
                    sfc = T(24.0) * (slj6 - T(2.0) * slj6 * slj6) / sdij_sq
                    f = sfc * vi1
                    forces[flane + i1] += f
                    forces[flane + j1] -= f
                end
            end)
        end
        # println("")

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
                flane4 = VecRange{4}(0)  # force lane
                flane  = VecRange{4}(0)  # force lane
                flanei = VecRange{4}(1)  # force lane
                forces = view(state.forces, :)
                idx = Vec{8,Int}(Tuple(1:8))
                lmask = Vec{8,Bool}((true,true,true,true,false,false,false,false))
                rmask = Vec{8,Bool}((false,false,false,false,true,true,true,true))
                # fi = zeros(T,4)
                # fi = view(zeros(T,8), :)
            end
        end)
        #endregion FORCE_SECTION
    
        m1 = (0, 4,  8, 12,  1,  5,  9, 13)
        m2 = (2, 6, 10, 14,  3,  7, 11, 15)
        mx = (0, 1,  2,  3,  8,  9, 10, 11)
        my = (4, 5,  6,  7, 12, 13, 14, 15)
        mj = (0,1,2,3,4,5,6,7)

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
            vi = vgather(coords, Vec{8,Int}((0,1,2,3,0,1,2,3)) + i1)

            while ptr+7 <= ptr_stop
                j = ptr
                ptr += 8
    
                # j = 4*(j-1) + 1 = 4j-4+1 = 4j-3 = j<<2 - 3
                
                j1 = vlist.list[j  ]<<2 - 3
                j2 = vlist.list[j+1]<<2 - 3
                j3 = vlist.list[j+2]<<2 - 3
                j4 = vlist.list[j+3]<<2 - 3
                j5 = vlist.list[j+4]<<2 - 3
                j6 = vlist.list[j+5]<<2 - 3
                j7 = vlist.list[j+6]<<2 - 3
                j8 = vlist.list[j+7]<<2 - 3

                vi1 = vload(Vec{4,T}, coords, j1) - vi
                vi2 = vload(Vec{4,T}, coords, j2) - vi
                vi3 = vload(Vec{4,T}, coords, j3) - vi
                vi4 = vload(Vec{4,T}, coords, j4) - vi
                vi5 = vload(Vec{4,T}, coords, j5) - vi
                vi6 = vload(Vec{4,T}, coords, j6) - vi
                vi7 = vload(Vec{4,T}, coords, j7) - vi
                vi8 = vload(Vec{4,T}, coords, j8) - vi


                
                
                # mj = (0,1,2,3,8,9,10,11)
                vi12 = shufflevector(vi1, vi2, Val{mj})
                vi34 = shufflevector(vi3, vi4, Val{mj})
                vi56 = shufflevector(vi5, vi6, Val{mj})
                vi78 = shufflevector(vi7, vi8, Val{mj})

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

                # vs12_m1 = shufflevector(vi1, vi2, Val{m1})    # xi1, xi2, yi1, yi2
                # vs34_m1 = shufflevector(vi3, vi4, Val{m1})    # xi3, xi4, yi3, yi4
                # vs12_m2 = shufflevector(vi1, vi2, Val{m2})    # zi1, zi2, ?i1, ?i2
                # vs34_m2 = shufflevector(vi3, vi4, Val{m2})    # zi3, zi4, ?i3, ?i4
                
                # xij = shufflevector(vs12_m1, vs34_m1, Val{mx})   # xi1, xi2, xi3, xi4
                # yij = shufflevector(vs12_m1, vs34_m1, Val{my})   # yi1, yi2, yi3, yi4
                # zij = shufflevector(vs12_m2, vs34_m2, Val{mx})   # zi1, zi2, zi3, zi4
    
                dij_sq = xij*xij + yij*yij + zij*zij
                mask = dij_sq <= cutsq
                !any(mask) && continue
    
                ϵm = ϵ[mlane, mask]
    
                lj2 = σ / dij_sq
                lj6 = lj2*lj2*lj2
                energy += sum( ϵm*(lj6*lj6 - lj6) )
                
                # println(typeof(energy), " ", typeof(sum( ϵm*(lj6*lj6 - lj6) )))
                # println(typeof(lj2))

                # region FORCE_SECTION
                $(if Forces === true
                    quote
                        fc = ϵm * T(24.0) * (lj6 - T(2.0) * lj6 * lj6) / dij_sq
                        #if mask[1]
                            # f = fc[1] * vi12
                            v = fc[1] * shufflevector(vi12, Val{(0,1,2,3)})
                            fi += v
                            # forces[flane+i1] += v
                            forces[flane+j1] -= v
                        #end
                        #if mask[2]
                            # f = fc[2] * vi12
                            v = fc[2] * shufflevector(vi12, Val{(4,5,6,7)})
                            fi += v
                            # forces[flane+i1] += v
                            forces[flane+j2] -= v
                        #end
                        #if mask[3]
                            # f = fc[3] * vi34
                            v = fc[3] * shufflevector(vi34, Val{(0,1,2,3)})
                            fi += v
                            # forces[flane+i1] += v
                            forces[flane+j3] -= v
                        #end
                        #if mask[4]
                            # f = fc[4] * vi34
                            v = fc[4] * shufflevector(vi34, Val{(4,5,6,7)})
                            fi += v
                            # forces[flane+i1] += v
                            forces[flane+j4] -= v
                        #end
                        #if mask[5]
                            # f = fc[5] * vi56
                            v = fc[5] * shufflevector(vi56, Val{(0,1,2,3)})
                            fi += v
                            # forces[flane+i1] += v
                            forces[flane+j5] -= v
                        #end
                        #if mask[6]
                            # f = fc[6] * vi56
                            v = fc[6] * shufflevector(vi56, Val{(4,5,6,7)})
                            fi += v
                            # forces[flane+i1] += v
                            forces[flane+j6] -= v
                        #end
                        #if mask[7]
                            # f = fc[7] * vi78
                            v = fc[7] * shufflevector(vi78, Val{(0,1,2,3)})
                            fi += v
                            # forces[flane+i1] += v
                            forces[flane+j7] -= v
                        #end
                        #if mask[8]
                            # f = fc[8] * vi78
                            v = fc[8] * shufflevector(vi78, Val{(4,5,6,7)})
                            fi += v
                            # forces[flane+i1] += v
                            forces[flane+j8] -= v
                        #end
                    end
                end)
                #endregion FORCE_SECTION
    
            end # for i
            
            
            #region FORCE_SECTION
            $(if Forces === true
                :(forces[flane + i1] += fi)
                # quote
                #     forces[i1  ] += fi[1]# + fi[5]
                #     forces[i1+1] += fi[2]# + fi[6]
                #     forces[i1+2] += fi[3]# + fi[7]
                # end
            end)
            #endregion FORCE_SECTION
            
            
            # do remaining pairs (atom-by-atom)
            vii = vload(Vec{4,T}, coords, i1)
            for j = ptr:ptr_stop
                #println("$i  $(vlist.list[j])")
                j1 = vlist.list[j]<<2 - 3
                vv1 = vload(Vec{4,T}, coords, j1)
                vi1 = vv1 - vii
            
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
                        # forces[flane4 + i1] += f
                        forces[flane4 + j1] -= f
                    end
                end)
            end
            #println("")
            
        end # for i
    
    end # quote
    end # generated function