using Base.Cartesian
import ..Calculators.eval!

@generated function eval!(state::State{T}, bonds::Vector{HarmonicPotential{Int,2,T}}, ::Type{Val{DoF}}) where {T<:AbstractFloat, DoF}
println("generating for $DoF, $T")
quote
    energy = $(T(0))
    coords = state.x
    forces = state.f
    @inbounds for bond in bonds
        a1 = bond.key[1]
        a2 = bond.key[2]
        
        @nexprs 3 u -> v12_u = coords[a2,u] - coords[a1,u]
        d12 = sqrt(ProtoSyn.@dot(u, v12_u, v12_u))
        dr = d12 - bond.x0

        energy += bond.kf * dr * dr

        #region FORCE_SECTION
        $(if DoF === true
            quote
                fc = bond.kf * dr / d12
                @nexprs 3 u -> forces[a1, u] += fc * v12_u
                @nexprs 3 u -> forces[a2, u] -= fc * v12_u
            end
        end)
        #endregion

    end
    state.e[HarmonicPotential{Int,2,T}] = energy/2
    state
end
end




@generated function eval!(state::State{T}, angles::Vector{HarmonicPotential{Int,3,T}}, ::Type{Val{DoF}}) where {T<:AbstractFloat, DoF}
quote
    energy = $(T(0))
    coords = state.x
    forces = state.f
    
    @inbounds for angle in angles
        a1,a2,a3 = angle.key

        @nexprs 3 u -> v21_u = coords[a1,u] - coords[a2,u]
        @nexprs 3 u -> v23_u = coords[a3,u] - coords[a2,u]

        # 1.0/(v21⋅v21)
        inv_d21sq  = 1.0 / @dot(u, v21_u, v21_u)
        
        # 1.0/(v23⋅v23)
        inv_d23sq  = 1.0 / @dot(u, v23_u, v23_u)
        
        # cos(θ) = v21⋅v23/(|v21|*|v23|)
        ctheta = @dot(u, v21_u, v23_u)

        inv_d21xd23 = sqrt(inv_d21sq * inv_d23sq)
        ctheta *= inv_d21xd23
        
        dtheta  = acos(ctheta) - angle.x0
        energy += angle.kf * dtheta * dtheta
        
        #region FORCE_SECTION
        $(if DoF === true
            quote
                # prevent division by zero
                fc = max(1.0 - ctheta*ctheta, 1e-8)
                fc = dtheta * type.kf / sqrt(fc)
                
                @nexprs 3 u -> begin
                    f1 = fc*(ctheta * v21_u*inv_d21sq - v23_u*inv_d21xd23)
                    f3 = fc*(ctheta * v23_u*inv_d23sq - v21_u*inv_d21xd23)
                    forces[a1, u] -= f1
                    forces[a2, u] += (f1 + f3)
                    forces[a3, u] -= f3
                end
            end
        end)
        #endregion
            
    end
    state.e[HarmonicPotential{Int,3,T}] = energy/2
    state
end
end

