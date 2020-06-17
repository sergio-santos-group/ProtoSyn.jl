@generated function eval!(state::State{T}, items::Vector{CoulombLennardJonesPotential{Int,2,T}}, pname::Symbol, ::Type{Val{DoF}}) where {T<:AbstractFloat, DoF}
quote
    eLJ = $(T(0))
    eCoulomb = $(T(0))
    fSI = $(T(1389.35458))  # kJ mol−1 Å e−2
    
    coords = state.x
    forces = state.f
    
    @inbounds for item in items
        a1 = item.key[1]
        a2 = item.key[2]
        
        @nexprs 3 u -> v12_u = coords[a2,u] - coords[a1,u]
        d12sq = @dot(u, v12_u, v12_u)
        
        lj2 = item.σ*item.σ/d12sq
        lj6 = lj2*lj2*lj2

        eLJ += item.ϵ * (lj6 * lj6 - lj6)

        ecoul = fSI * item.q / sqrt(d12sq)
        eCoulomb += ecoul

        #region FORCE_SECTION
        $(if DoF === true
            quote
                fc = (24 * item.ϵ * (lj6 - 2 * lj6 * lj6) - ecoul) / d12sq
                @nexprs 3 u -> forces[a1, u] += fc * v12_u
                @nexprs 3 u -> forces[a2, u] -= fc * v12_u
            end
        end)
        #endregion

    end
    #state.e[:CoulombPair] = eCoulomb
    #state.e[:LJPair] = eLJ
    
    
    # state.e[:($(name(eltype(items), "LJ")))] = eLJ
    # state.e[:($(name(eltype(items), "Coul")))] = eCoulomb

    setenergy!(state, :LJPair, eLJ)
    setenergy!(state, :CoulombPair, eCoulomb)
    state
end
end

@generated function eval!(state::State{T}, bonds::Vector{CoulombLennardJonesPotential{Int,1,T}}, pname::Symbol, ::Type{Val{DoF}}) where {T<:AbstractFloat, DoF}
quote
    #state.e[$(println("ola"); :(name(eltype(bonds))))] = $(T(0))
    # $(T(0))
    state
end
end