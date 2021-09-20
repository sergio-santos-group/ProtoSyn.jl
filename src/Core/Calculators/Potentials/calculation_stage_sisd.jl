# Main function
@generated function resolve_calculation(::Type{ProtoSyn.SISD_0},
    pose::Pose,
    potential::Function,
    update_forces::Bool,
    verlet_list::Nothing,
    coords::Vector{T},
    mask::MaskMap) where {T <: AbstractFloat}

    quote
        n_atoms = length(coords) ÷ 3
        forces   = zeros(T, n_atoms, n_atoms, 3)
        energies = zeros(T, n_atoms, n_atoms)
        @inbounds for i in 1:n_atoms - 1
            i1 = i<<2 - (2 + i) # Conversion from atom number to array position
            @nexprs 3 u -> ri_u = coords[i1 + u - 1]
            
            @inbounds for j = i+1:n_atoms
                j1 = j<<2 - (2 + j) # Conversion from atom number to array position
                @nexprs 3 u -> rij_u = coords[j1 + u - 1] - ri_u
                dij = sqrt(@reduce 3 (+) u -> rij_u * rij_u)
                @nexprs 3 u -> rij_u = rij_u / dij # normalization
                energies[j, i], (forces[j, i, 1], forces[j, i, 2], forces[j, i, 3]), (forces[i, j, 1], forces[i, j, 2], forces[i, j, 3]) = potential(dij, v = (rij_1, rij_2, rij_3))
            end
        end
        
        return resolve_mask(ProtoSyn.SISD_0, pose, energies, forces, mask)
    end
end

# TODO: Verlet list support