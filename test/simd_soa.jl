@generated function simd_soa(state::State{T}, cutoff::Float64, ::Type{Val{Forces}}, ::Type{Vec{N,T}}) where {Forces,N,T<:AbstractFloat}
    
    @assert isa(Forces, Bool) "Forces must be a boolean"

    fbody = quote
        
        σ = T(1)
        ϵ = ones(T, N)
        energy = T(0)
        natoms = state.size
        mlane = VecRange{N}(1)  # mask lane
        cutsq = convert(T, cutoff*cutoff)
        vcutsq = Vec{N,T}(cutoff*cutoff)

        @nexprs 3 u -> x_u = view(state.coords, :, u)
        
        #region FORCE_SECTION
        $(if Forces === true
            quote
                fill!(state.forces, T(0))
                @nexprs 3 u -> f_u = view(state.forces, :, u)
                flane = VecRange{N}(0)  # force lane
            end
        end)
        #endregion FORCE_SECTION
        
        @inbounds for i = 1:natoms-1
            
            @nexprs 3 u -> ri_u = Vec{N,T}(x_u[i])
            
            #region FORCE_SECTION
            $(if Forces === true
                quote
                    @nexprs 3 u ->  fi_u = Vec{N,T}(0)
                    @nexprs 3 u -> sfi_u = T(0)
                end
            end)
            #endregion FORCE_SECTION
            
            jj = i + 1
            while jj+N-1 <= natoms
                j = jj
                jj += N
                
                @nexprs 3 u -> rij_u = vload(Vec{N,T}, x_u, j) - ri_u
                dij_sq = @reduce 3 (+) u -> rij_u*rij_u
                
                mask = dij_sq <= vcutsq
                !any(mask) && continue
                
                ϵm = ϵ[mlane, mask]

                lj2 = σ / dij_sq
                lj6 = lj2*lj2*lj2
                energy += sum( ϵm*(lj6*lj6 - lj6) )
                #energy += sum( (lj6*lj6 - lj6) )

                #region FORCE_SECTION
                $(if Forces === true
                    quote
                        fc = ϵm * T(24.0) * (lj6 - T(2.0) * lj6 * lj6) / dij_sq
                        #fc = T(24.0) * (lj6 - T(2.0) * lj6 * lj6) / dij_sq
                        lanej = flane + j
                        @nexprs 3 u -> begin
                            ftmp = fc * rij_u
                            f_u[lanej, mask] -= ftmp
                            fi_u += ftmp
                        end
                    end
                end)
                #endregion FORCE_SECTION
            end

            # do remaining pairs
            @nexprs 3 u -> sri_u = x_u[i]
            for j=jj:natoms

                @nexprs 3 u -> srij_u = x_u[j] - sri_u
                sdij_sq = @reduce 3 (+) u -> srij_u * srij_u
                (sdij_sq > cutsq) && continue

                slj2 = σ / sdij_sq
                slj6 = slj2*slj2*slj2
                energy += sum( slj6*slj6 - slj6 )
                
                #region FORCE_SECTION
                $(if Forces === true
                    quote
                        sfc = T(24.0) * (slj6 - T(2.0) * slj6 * slj6) / sdij_sq
                        # @nexprs 3 u -> f_u[i] += sfc * srij_u
                        @nexprs 3 u -> sfi_u += sfc * srij_u
                        @nexprs 3 u -> f_u[j] -= sfc * srij_u
                    end
                end)
                #endregion FORCE_SECTION
            end

            #region FORCE_SECTION
            $(if Forces === true
               :(@nexprs 3 u -> f_u[i] += sum(fi_u) + sfi_u)
            end)
           #endregion FORCE_SECTION
            
        end
        energy
        #lm
    end
    fbody
end