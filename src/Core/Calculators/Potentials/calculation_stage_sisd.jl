# Main function
@generated function resolve_calculation(::Type{ProtoSyn.SISD_0},
    pose::Pose,
    potential::Function,
    update_forces::Bool,
    verlet_list::Nothing,
    coords::Vector{T},
    mask::MaskMap,
    indexes::Vector{Int}) where {T <: AbstractFloat}

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

                if !update_forces
                    energies[j, i] = potential(dij, qi = pose.state[indexes[i]].δ, qj = pose.state[indexes[j]].δ)
                else
                    @nexprs 3 u -> rij_u = rij_u / dij # normalization
                    energies[j, i], (forces[j, i, 1], forces[j, i, 2], forces[j, i, 3]), (forces[i, j, 1], forces[i, j, 2], forces[i, j, 3]) = potential(dij, v = (rij_1, rij_2, rij_3), qi = pose.state[i].δ, qj = pose.state[j].δ)
                end
            end
        end
        
        return resolve_mask(ProtoSyn.SISD_0, pose, energies, forces, update_forces, mask)
    end
end

@generated function resolve_calculation(::Type{ProtoSyn.SISD_0},
    pose::Pose,
    potential::Function,
    update_forces::Bool,
    verlet_list::VerletList,
    coords::Vector{T},
    mask::MaskMap,
    indexes::Vector{Int}) where {T <: AbstractFloat}

    quote

        n_atoms = length(coords) ÷ 3
        verlet_list.size !== n_atoms && @warn "Verlet list size and number of selected atoms do not match"
        forces   = zeros(T, n_atoms, n_atoms, 3)
        energies = zeros(T, n_atoms, n_atoms)

        for i in 1:n_atoms - 1
            ptr = verlet_list.offset[i]
            if verlet_list.list[ptr] < 1
                continue
            end
            ptr_stop = verlet_list.offset[i+1]-2

            i1 = i<<2 - (2 + i) # Conversion from atom number to array position
            @nexprs 3 u -> ri_u = coords[i1 + u - 1]
            
            while ptr <= ptr_stop
                j = verlet_list.list[ptr]

                j1 = j<<2 - (2 + j) # Conversion from atom number to array position
                @nexprs 3 u -> rij_u = coords[j1 + u - 1] - ri_u
                dij = sqrt(@reduce 3 (+) u -> rij_u * rij_u)

                if !update_forces
                    energies[j, i] = potential(dij)
                else
                    @nexprs 3 u -> rij_u = rij_u / dij # normalization
                    energies[j, i], (forces[j, i, 1], forces[j, i, 2], forces[j, i, 3]), (forces[i, j, 1], forces[i, j, 2], forces[i, j, 3]) = potential(dij, v = (rij_1, rij_2, rij_3), qi = pose.state[indexes[i]].δ, qj = pose.state[indexes[j]].δ)
                end

                ptr += 1
            end
        end
        
        return resolve_mask(ProtoSyn.SISD_0, pose, energies, forces, update_forces, mask)
    end
end