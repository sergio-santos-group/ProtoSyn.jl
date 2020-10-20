function simd(state::State{T}) where {T <: AbstractFloat}

    natoms          = state.size
    distance_matrix = zeros(natoms, natoms)
    remaining_mask  = Vec{4, T}((1, 1, 1, 0))

    # For the coords we need to add an extra 0, because the last coord will
    # be loaded in a 4-wide vector for SIMD calculation. The last value is
    # then ignored.
    coords = push!(state.x[:], T(0))

    @inbounds for i = 1:natoms-1

        # Atom number to atom position conversion
        _i = i << 2 - (2 + i)
        vi = vload(Vec{4, T}, coords, _i) # Load XYZ (consecutive)

        @inbounds for j = (i+1):natoms
            
            # Atom number to atom position conversion
            _j = j << 2 - (2 + j)

            # Calculate distance, while ignoring the last value in the
            # 4-wide SIMD vector
            rij = (vload(Vec{4, T}, coords, _j) - vi) * remaining_mask   # xi1, yi1, zi1, ?i1
            distance_matrix[i, j] = sqrt(sum(rij*rij))
        end # for j
    end # for i

    return distance_matrix
end


function simd(state::State{T}, verlet_list::VerletList) where {T <: AbstractFloat}

    natoms          = state.size
    distance_matrix = zeros(natoms, natoms)
    remaining_mask  = Vec{4, T}((1, 1, 1, 0))

    # For the coords we need to add an extra 0, because the last coord will
    # be loaded in a 4-wide vector for SIMD calculation. The last value is
    # then ignored.
    coords = push!(state.x[:], T(0))

    @inbounds for i = 1:natoms-1

        ptr = verlet_list.offset[i]
        if verlet_list.list[ptr] < 1
            continue
        end
        ptr_stop = verlet_list.offset[i+1]-2

        # Atom number to atom position conversion
        i1 = (i<<1) + (i-2)
        vi = vload(Vec{4, T}, coords, i1) # Load XYZ (consecutive)

        while ptr <= ptr_stop
            j = ptr
            ptr += 1
            
            _j = verlet_list.list[j]<<2 - (2 + verlet_list.list[j])
            
            # Calculate distance, while ignoring the last value in the
            # 4-wide SIMD vector
            rij = (vload(Vec{4, T}, coords, _j) - vi) * remaining_mask # xi1, yi1, zi1, wi1
            
            distance_matrix[i, verlet_list.list[j]] = sqrt(sum(rij * rij))
        end # for j
    end # for i

    return distance_matrix
end