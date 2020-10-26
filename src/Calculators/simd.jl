export distance_matrix_simd

"""
    Calculators.distance_matrix_simd(coords::Vector{T}) where {T <: AbstractFloat}
    Calculators.distance_matrix_simd(state::State{T}) where {T <: AbstractFloat}
    Calculators.distance_matrix_simd(pose::Pose)

Return a distance matrix with the distance of all pairs of coordinates in
`coords` (this should be a Matrix{T} in AoS format), in the CPU using SIMD
architecture. Instead of a Vector{T} `coords`, a `State` or `Pose` can be
provided, in which case the coordinates considered are all the existent in the
State or Pose.state, respectively.

    Calculators.distance_matrix_simd(pose::Pose, selection::ProtoSyn.AbstractSelection)

Return a distance matrix with the distance of all pairs of atoms in the Pose
`pose` who are included in the given `selection`, in the CPU using SIMD
architecture.

    Calculators.distance_matrix_simd(coords::Vector{T}, vlist::VerletList) where {T <: AbstractFloat}
    Calculators.distance_matrix_simd(state::State{T}, vlist::VerletList) where {T <: AbstractFloat}
    Calculators.distance_matrix_simd(pose::Pose, vlist::VerletList)

Return a distance matrix with the distance of all pairs of coordinates in the
`VerletList` `vlist`, in the CPU using SIMD architecture. Check `VerletList` for 
a more in-depth look at how Verlet lists work. Instead of a Vector{T} `coords`,
a `State` or `Pose` can be provided, in which case the coordinates considered
are existent in the State or Pose.state, respectively. _Note:_ Selections can
still be applied when using Verlet lists, but need to be applied when updating
the lists themselves. Check `VerletList` for a more in-depth look at how Verlet
lists work.

# See also

`distance_matrix_serial`, `distance_matrix_cuda`

# Examples
```jldoctest
julia> Calculators.distance_matrix_simd(pose.state.x)
N×N CUDA.CuArray{Float64,2}:
...

julia> Calculators.distance_matrix_simd(pose, an"CA")
N×N CUDA.CuArray{Float64,2}:
...

julia> Calculators.distance_matrix_simd(pose, verlet_list)
N×N CUDA.CuArray{Float64,2}:
...
```
"""
function distance_matrix_simd(coords::Vector{T}) where {T <: AbstractFloat}
    # coords must be in AoS format

    natoms          = trunc(Int64, length(coords)/3)
    distance_matrix = zeros(natoms, natoms)
    remaining_mask  = Vec{4, T}((1, 1, 1, 0))

    # For the coords we need to add an extra 0, because the last coord will
    # be loaded in a 4-wide vector for SIMD calculation. The last value is
    # then ignored.
    _coords = push!(copy(coords), T(0))

    @inbounds for i = 1:natoms-1

        # Atom number to atom position conversion
        _i = i << 2 - (2 + i)
        vi = vload(Vec{4, T}, _coords, _i) # Load XYZ (consecutive)

        @inbounds for j = (i+1):natoms
            
            # Atom number to atom position conversion
            _j = j << 2 - (2 + j)

            # Calculate distance, while ignoring the last value in the
            # 4-wide SIMD vector
            rij = (vload(Vec{4, T}, _coords, _j) - vi) * remaining_mask   # xi1, yi1, zi1, ?i1
            distance_matrix[i, j] = sqrt(sum(rij*rij))
        end # for j
    end # for i

    return distance_matrix
end

distance_matrix_simd(pose::Pose) = distance_matrix_simd(pose.state.x[:])
distance_matrix_simd(state::State{T}) where {T <: AbstractFloat} = begin
    distance_matrix_simd(state.x[:])
end


distance_matrix_simd(pose::Pose, selection::ProtoSyn.AbstractSelection) = begin

    s = selection(pose, gather = true)
    coords = zeros(eltype(pose.state), 3, length(s))
    for (i, atom) in enumerate(s)
        coords[:, i] = pose.state[atom].t
    end

    return distance_matrix_simd(coords[:])
end


function distance_matrix_simd(coords::Vector{T}, verlet_list::VerletList) where {T <: AbstractFloat}

    natoms          = trunc(Int64, length(coords)/3)
    distance_matrix = zeros(natoms, natoms)
    remaining_mask  = Vec{4, T}((1, 1, 1, 0))

    # For the coords we need to add an extra 0, because the last coord will
    # be loaded in a 4-wide vector for SIMD calculation. The last value is
    # then ignored.
    _coords = push!(copy(coords), T(0))

    @inbounds for i = 1:natoms-1

        ptr = verlet_list.offset[i]
        if verlet_list.list[ptr] < 1
            continue
        end
        ptr_stop = verlet_list.offset[i + 1] - 2

        # Atom number to atom position conversion
        i1 = (i<<1) + (i-2)
        vi = vload(Vec{4, T}, _coords, i1) # Load XYZ (consecutive)

        while ptr <= ptr_stop
            j = ptr
            ptr += 1
            
            _j = verlet_list.list[j]<<2 - (2 + verlet_list.list[j])
            
            # Calculate distance, while ignoring the last value in the
            # 4-wide SIMD vector
            rij = (vload(Vec{4, T}, _coords, _j) - vi) * remaining_mask # xi1, yi1, zi1, wi1
            
            distance_matrix[i, verlet_list.list[j]] = sqrt(sum(rij * rij))
        end # for j
    end # for i

    return distance_matrix
end

distance_matrix_simd(state::State{T}, vlist::VerletList) where {T <: AbstractFloat} = begin
    distance_matrix_simd(state.x[:], vlist)
end

distance_matrix_simd(pose::Pose, vlist::VerletList) = begin
    distance_matrix_simd(pose.state.x[:], vlist)
end