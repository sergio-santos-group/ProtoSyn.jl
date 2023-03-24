using SIMD

# ! Note:
# When applying SIMD optimizations, it is possible to use 8 Float32 numbers in a
# single modern register (256 bits). With the right vector suffle operations
# this allows for the quick calculation of 8 inter-atomic euclidean distances in
# a single vector. However, SIMD vectors are not iterable → a single operation
# must be applied to all the data. For this reason, conditional potentials, such
# as the flat-bottom potential, are hard to apply in SIMD mode (the actual
# potential to be applied depends on the distance). An option, in these cases,
# is to use a mask for each of the domains in the flat-bottom potential,
# applying each sub-potential only to the distances in that given domain (the
# un-masked distances). In the current organization of ProtoSyn, that would
# require the existance of multiple potential methods → one whose argument is a
# single float (as in SISD_0 and CUDA_2 modes) and another whose argument is a
# SIMD vector. However, when looping over the atom interactions, in most use
# cases the number of interactions with any given atom is not a multiple of 8,
# and there there is a number of "left over" interaction to calculate. The most
# simple approach is to deal with these interactions one-at-a-time. Using this
# approach instead of employing the SIMD potential method, we now require the
# original 1-on-1 SIMD/CUDA potential. Providing both potentials to the
# `resolve_calculation` bellow is not feasible in the current architecture of
# ProtoSyn. Another option is to include "fake interactions" to fill the
# remaining slots in an 8 float32 register, which are later ignored. Moreover,
# for a correct calculation of the forces felt on each interacting atom, the
# displacement vector v must be normalized. Since the vector length is being
# calculated in the calculation stage, the potential function already receives a
# normalized vector. For modularity, non pre-normalized vectors can be used,
# with the addition of the `is_normalized` flag. In the case of it being set to
# false, the normalization is done during the potential application. 

@generated function resolve_calculation(::Type{ProtoSyn.SIMD_1},
    pose::Pose,
    potential::Function,
    update_forces::Bool,
    verlet_list::Nothing,
    coords::Vector{T},
    mask::MaskMap,
    indexes::Vector{Int}) where {T <: AbstractFloat}

    quote
        coords   = vcat(pose.state.x.coords[:], [0.0, 0.0, 0.0])
        n_atoms  = length(indexes)
        forces   = zeros(T, n_atoms, n_atoms, 3)
        energies = zeros(T, n_atoms, n_atoms)
        _mask = Vec{4, T}((1, 1, 1, 0))

        if update_forces
            arr1 = [T(0.0), T(0.0), T(0.0), T(0.0)] # Temporary save
            arr2 = [T(0.0), T(0.0), T(0.0), T(0.0)] # Temporary save
        end

        @inbounds @simd for i in 1:n_atoms - 1
            _i = indexes[i]
            i1 = _i<<2 - (2 + _i) # Conversion from atom number to array position
            vi = vload(Vec{4, T}, coords, i1)
            @inbounds @simd for j = i+1:n_atoms
                _j = indexes[j]
                j1 = _j<<2 - (2 + _j) # Conversion from atom number to array position
                rij = (vload(Vec{4, T}, coords, j1) - vi) * _mask
                dij = sqrt(sum(rij * rij))

                if !update_forces
                    energies[j, i] = potential(dij, qi = pose.state[indexes[i]].δ, qj = pose.state[indexes[j]].δ)
                else
                    rij = rij / dij
                    energies[j, i], f1, f2 = potential(dij, v = rij, qi = pose.state[indexes[i]].δ, qj = pose.state[indexes[j]].δ)
                    vstore(f1, arr1, 1)
                    vstore(f2, arr2, 1)
                    forces[j, i, 1] = f1[1]
                    forces[j, i, 2] = f1[2]
                    forces[j, i, 3] = f1[3]
                    forces[i, j, 1] = f2[1]
                    forces[i, j, 2] = f2[2]
                    forces[i, j, 3] = f2[3]
                end
            end
        end
        
        return resolve_mask(ProtoSyn.SIMD_1, pose, energies, forces, update_forces, mask)
    end
end


# With VerletList
@generated function resolve_calculation(::Type{ProtoSyn.SIMD_1},
    pose::Pose,
    potential::Function,
    update_forces::Bool,
    verlet_list::VerletList,
    coords::Vector{T},
    mask::MaskMap,
    indexes::Vector{Int}) where {T <: AbstractFloat}

    quote
        n_atoms  = length(coords) ÷ 3
        forces   = zeros(T, n_atoms, n_atoms, 3)
        energies = zeros(T, n_atoms, n_atoms)
        arr1 = [T(0.0), T(0.0), T(0.0)] # Temporary save
        arr2 = [T(0.0), T(0.0), T(0.0)] # Temporary save
        _mask = Vec{4, T}((1, 1, 1, 0))

        @inbounds for i in 1:n_atoms - 1

            ptr = verlet_list.offset[i]
            if verlet_list.list[ptr] < 1
                continue
            end
            ptr_stop = verlet_list.offset[i+1]-2

            i1 = i<<2 - (2 + i) # Conversion from atom number to array position
            vi = vload(Vec{4, T}, coords, i1)
            
            while ptr <= ptr_stop
                j = verlet_list.list[ptr]
                j1 = j<<2 - (2 + j) # Conversion from atom number to array position
                rij = (vload(Vec{4, T}, coords, j1) - vi) * _mask
                dij = sqrt(sum(rij * rij))

                if !update_forces
                    energies[j, i] = potential(dij)
                else
                    rij = rij / dij # normalization
                    energies[j, i], f1, f2 = potential(dij, v = rij, qi = pose.state[indexes[i]].δ, qj = pose.state[indexes[j]].δ)
                    vstore(f1, arr1, 1)
                    vstore(f2, arr2, 1)
                    forces[j, i, 1] = f1[1]
                    forces[j, i, 2] = f1[2]
                    forces[j, i, 3] = f1[3]
                    forces[i, j, 1] = f2[1]
                    forces[i, j, 2] = f2[2]
                    forces[i, j, 3] = f2[3]
                end

                ptr += 1
            end
        end
        
        return resolve_mask(ProtoSyn.SIMD_1, pose, energies, forces, update_forces, mask)
    end
end