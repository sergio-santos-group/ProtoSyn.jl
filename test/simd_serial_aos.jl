@generated function serial_aos(state::State{T}, cutoff::Float64, ::Type{Val{DOF}}) where {DOF, T <: AbstractFloat}
    # println("generating serial3N with DOF = $DOF")
    fbody = quote
        coords = state.coords
        forces = state.forces
        natoms = state.size

        cutsq = convert(T, cutoff*cutoff)   # squared cutoff
        energy = T(0)                       # total energy
        σ = T(1)                            # LJ sigma
        
        @inbounds for i = 1:natoms-1

            # load coordinates for the i-th atom
            @nexprs 3 u -> ri_u = coords[u,i]
            
            #region FORCE_SECTION <><><><><><><><><><><><><><><><><><><><><><><>
            # zero the force accummulator for the i-th atom
            $(if DOF === true
                :(@nexprs 3 u -> fi_u = zero(T))
            end)
            #endregion  <><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            
            for j = i+1:natoms

                # load coordinates for the j-th atom
                # and calculate the ij vector
                @nexprs 3 u -> rij_u = coords[u, j] - ri_u
                
                # calculate the squared distance. Skip
                # if greater than cutoff
                dij_sq = @reduce 3 (+) u -> rij_u*rij_u
                (dij_sq > cutsq) && continue
                

                # LJ potential
                lj2 = σ/dij_sq
                lj6 = lj2*lj2*lj2
                energy += lj6*lj6 - lj6
                
                #region FORCE_SECTION <><><><><><><><><><><><><><><><><><><><>
                $(if DOF === true
                    quote
                        fc = T(24) * (lj6 - T(2) * lj6 * lj6) / dij_sq
                        @nexprs 3 u -> begin
                            fi_u += fc * rij_u
                            forces[u, j] -= fc * rij_u
                        end
                    end
                end)
                #endregion  <><><><><><><><><><><><><><><><><><><><><><><><><><>
            
            end
            #region FORCE_SECTION <><><><><><><><><><><><><><><><><><><><><><><>
            # update forces for the i-th atom
            $(if DOF === true
                :(@nexprs 3 u -> forces[u, i] += fi_u)
            end)
            #endregion  <><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        end
        energy
    end
    fbody
end







@generated function serial_aos(state::State{T}, vlist::VerletList, ::Type{Val{DOF}}) where {DOF, T <: AbstractFloat}
    # println("generating serial3N with DOF = $DOF")
    fbody = quote
        coords = state.coords
        forces = state.forces
        natoms = state.size

        cutsq = convert(T, vlist.cutoff*vlist.cutoff)   # squared cutoff
        energy = T(0)                       # total energy
        σ = T(1)                            # LJ sigma
        
        @inbounds for i = 1:natoms-1

            # load coordinates for the i-th atom
            @nexprs 3 u -> ri_u = coords[u,i]
            
            #region FORCE_SECTION <><><><><><><><><><><><><><><><><><><><><><><>
            # zero the force accummulator for the i-th atom
            $(if DOF === true
                :(@nexprs 3 u -> fi_u = zero(T))
            end)
            #endregion  <><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            
            ptr = vlist.offset[i]
            while vlist.list[ptr] > 0
                j = vlist.list[ptr]
                ptr += 1

                # load coordinates for the j-th atom
                # and calculate the ij vector
                @nexprs 3 u -> rij_u = coords[u, j] - ri_u
                
                # calculate the squared distance. Skip
                # if greater than cutoff
                dij_sq = @reduce 3 (+) u -> rij_u*rij_u
                (dij_sq > cutsq) && continue
                

                # LJ potential
                lj2 = σ/dij_sq
                lj6 = lj2*lj2*lj2
                energy += lj6*lj6 - lj6
                
                #region FORCE_SECTION <><><><><><><><><><><><><><><><><><><><>
                $(if DOF === true
                    quote
                        fc = T(24) * (lj6 - T(2) * lj6 * lj6) / dij_sq
                        @nexprs 3 u -> begin
                            fi_u += fc * rij_u
                            forces[u, j] -= fc * rij_u
                        end
                    end
                end)
                #endregion  <><><><><><><><><><><><><><><><><><><><><><><><><><>
            
            end
            #region FORCE_SECTION <><><><><><><><><><><><><><><><><><><><><><><>
            # update forces for the i-th atom
            $(if DOF === true
                :(@nexprs 3 u -> forces[u, i] += fi_u)
            end)
            #endregion  <><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        end
        energy
    end
    fbody
end