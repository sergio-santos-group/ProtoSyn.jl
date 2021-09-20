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

# Main function
# function resolve_calculation(::Type{ProtoSyn.SIMD_1},
#     pose::Pose,
#     potential::Function,
#     update_forces::Bool,
#     verlet_list::Nothing,
#     coords::Vector{T},
#     mask::MaskMap) where {T <: AbstractFloat}

#     n_atoms = length(coords) ÷ 3
#     coords  = vcat(coords, zeros(T, 24)) # fake interactions

#     # ! This is using N (register width) = 8
#     coords_size = length(coords) ÷ 3
#     forces   = zeros(T, coords_size * coords_size * 3)
#     energies = zeros(T, coords_size * coords_size)

#     mlane = VecRange{8}(1)
#     remaining_mask = Vec{4, T}((1, 1, 1, 0))

#     m1 = (0, 4,  8, 12,  1,  5,  9, 13)
#     m2 = (2, 6, 10, 14,  3,  7, 11, 15)
#     mx = (0, 1,  2,  3,  8,  9, 10, 11)
#     my = (4, 5,  6,  7, 12, 13, 14, 15)
#     mj = (0, 1,  2,  3,  4,  5,  6,  7)

#     @inbounds for i in 1:(n_atoms - 1)
#         i1 = i<<2 - (2 + i) # Conversion from atom number to array position
#         vi = vload(Vec{4, T}, coords, i1)
#         j  = i
#         while j < n_atoms
#             @nexprs 8 u ->  j_u = ((j+u)<<1) + (j+u-2) # Atom number to coordinate position
#             @nexprs 8 u -> vi_u = vload(Vec{4,T}, coords, j_u) - vi    # xi1, yi1, zi1, wi1
#             vi12 = shufflevector(vi_1, vi_2, Val(mj))           # xi1, yi1, zi1, wi1, xi2, yi2, zi2, wi2
#             vi34 = shufflevector(vi_3, vi_4, Val(mj))           # xi3, yi3, zi3, wi3, xi4, yi4, zi4, wi4
#             vi56 = shufflevector(vi_5, vi_6, Val(mj))           # xi5, yi5, zi5, wi5, xi6, yi6, zi6, wi6
#             vi78 = shufflevector(vi_7, vi_8, Val(mj))           # xi7, yi7, zi7, wi7, xi8, yi8, zi8, wi8
        
#             vs1234_m1 = shufflevector(vi12, vi34, Val(m1))      # xi1, xi2, xi3, xi4, yi1, yi2, yi3, yi4
#             vs5678_m1 = shufflevector(vi56, vi78, Val(m1))      # xi5, xi6, xi7, xi8, yi5, yi6, yi7, yi8
#             vs1234_m2 = shufflevector(vi12, vi34, Val(m2))      # zi1, zi2, zi3, zi4, wi1, wi2, wi3, wi4
#             vs5678_m2 = shufflevector(vi56, vi78, Val(m2))      # zi5, zi6, zi7, zi8, wi5, wi6, wi7, wi8

#             xij = shufflevector(vs1234_m1, vs5678_m1, Val(mx))  # xi1, xi2, xi3, xi4, xi5, xi6, xi7, xi8
#             yij = shufflevector(vs1234_m1, vs5678_m1, Val(my))  # yi1, yi2, yi3, yi4, yi5, yi6, yi7, yi8
#             zij = shufflevector(vs1234_m2, vs5678_m2, Val(mx))  # zi1, zi2, zi3, zi4, zi5, zi6, zi7, zi8

#             dij = sqrt(xij*xij + yij*yij + zij*zij)
#             xij /= dij    # normalize
#             yij /= dij    # normalize
#             zij /= dij    # normalize
#             e, f1, f2 = potential(dij, v = [xij, yij, zij])

#             # Store energy
#             vstore(e, energies, (i-1) * n_atoms + j + 1)

#             # Store forces
#             vstore(f1[1], forces, (i-1) * coords_size + j + 1)
#             vstore(f1[2], forces, (i-1) * coords_size + j + 1 + (coords_size * coords_size))
#             vstore(f1[3], forces, (i-1) * coords_size + j + 1 + (coords_size * coords_size * 2))

#             # Automating idx definition is too slow: manual attribution is faster
#             idx1 = Vec((j * coords_size + i, (j + 1) * coords_size + i, (j + 2) * coords_size + i, (j + 3) * coords_size + i, (j + 4) * coords_size + i, (j + 5) * coords_size + i, (j + 6) * coords_size + i, (j + 7) * coords_size + i))
#             idx2 = Vec((j * coords_size + i + (coords_size * coords_size), (j + 1) * coords_size + i + (coords_size * coords_size), (j + 2) * coords_size + i + (coords_size * coords_size), (j + 3) * coords_size + i + (coords_size * coords_size), (j + 4) * coords_size + i + (coords_size * coords_size), (j + 5) * coords_size + i + (coords_size * coords_size), (j + 6) * coords_size + i + (coords_size * coords_size), (j + 7) * coords_size + i + (coords_size * coords_size)))
#             idx3 = Vec((j * coords_size + i + (coords_size * coords_size * 2), (j + 1) * coords_size + i + (coords_size * coords_size * 2), (j + 2) * coords_size + i + (coords_size * coords_size * 2), (j + 3) * coords_size + i + (coords_size * coords_size * 2), (j + 4) * coords_size + i + (coords_size * coords_size * 2), (j + 5) * coords_size + i + (coords_size * coords_size * 2), (j + 6) * coords_size + i + (coords_size * coords_size * 2), (j + 7) * coords_size + i + (coords_size * coords_size * 2)))
#             vscatter(f2[1], forces, idx1)
#             vscatter(f2[2], forces, idx2)
#             vscatter(f2[3], forces, idx3)

#             if j + 8 <= n_atoms
#                 j += 8 # Do the next N atoms
#             else
#                 break
#             end
#         end
#     end
    
#     energies = reshape(energies[1:(n_atoms * n_atoms)], (n_atoms, n_atoms))
#     forces   = reshape(forces, (coords_size, :))[1:n_atoms, :]
#     forces   = reshape(forces, (n_atoms, coords_size, :))[:, 1:n_atoms, :]
#     return energies, forces
# end

@generated function resolve_calculation(::Type{ProtoSyn.SIMD_1},
    pose::Pose,
    potential::Function,
    update_forces::Bool,
    verlet_list::Nothing,
    coords::Vector{T},
    mask::MaskMap) where {T <: AbstractFloat}

    quote
        n_atoms  = length(coords) ÷ 3
        forces   = zeros(T, n_atoms, n_atoms, 3)
        energies = zeros(T, n_atoms, n_atoms)
        arr1 = [T(0.0), T(0.0), T(0.0)] # Temporary save
        arr2 = [T(0.0), T(0.0), T(0.0)] # Temporary save
        _mask = Vec{4, T}((1, 1, 1, 0))

        @inbounds @simd for i in 1:n_atoms - 1
            i1 = i<<2 - (2 + i) # Conversion from atom number to array position
            vi = vload(Vec{4, T}, coords, i1)
            
            @inbounds @simd for j = i+1:n_atoms
                j1 = j<<2 - (2 + j) # Conversion from atom number to array position
                rij = (vload(Vec{4, T}, coords, j1) - vi) * _mask
                dij = sqrt(sum(rij * rij))
                rij = rij / dij # normalization
                energies[j, i], f1, f2 = potential(dij, v = rij)
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
        
        return resolve_mask(ProtoSyn.SIMD_1, pose, energies, forces, mask)
    end
end

# TODO: Verlet list support