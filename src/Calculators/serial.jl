export distance_matrix_serial

"""
    Calculators.distance_matrix_serial(coords::Matrix{T}) where {T <: AbstractFloat}
    Calculators.distance_matrix_serial(state::State{T}) where {T <: AbstractFloat}
    Calculators.distance_matrix_serial(pose::Pose)

Return a distance matrix with the distance of all pairs of coordinates in
`coords` (this should be a Matrix{T} in AoS format), in the CPU using SISD
architecture. Instead of a Matrix{T} `coords`, a `State` or `Pose` can be
provided, in which case the coordinates considered are all the existent in the
State or Pose.state, respectively.

    Calculators.distance_matrix_serial(pose::Pose, selection::ProtoSyn.AbstractSelection)

Return a distance matrix with the distance of all pairs of atoms in the Pose
`pose` who are included in the given `selection`, in the CPU using SISD
architecture.

    Calculators.distance_matrix_serial(coords::Matrix{T}, vlist::VerletList) where {T <: AbstractFloat}
    Calculators.distance_matrix_serial(state::State{T}, vlist::VerletList) where {T <: AbstractFloat}
    Calculators.distance_matrix_serial(pose::Pose, vlist::VerletList)

Return a distance matrix with the distance of all pairs of coordinates in the
`VerletList` `vlist`, in the CPU using SISD architecture. Check `VerletList` for 
a more in-depth look at how Verlet lists work. Instead of a Matrix{T} `coords`,
a `State` or `Pose` can be provided, in which case the coordinates considered
are existent in the State or Pose.state, respectively. _Note:_ Selections can
still be applied when using Verlet lists, but need to be applied when updating
the lists themselves. Check `VerletList` for a more in-depth look at how Verlet
lists work.

# See also

`distance_matrix_simd`, `distance_matrix_cuda`

# Examples
```jldoctest
julia> Calculators.distance_matrix_serial(pose.state.x)
N×N CUDA.CuArray{Float64,2}:
    ...

julia> Calculators.distance_matrix_serial(pose, an"CA")
N×N CUDA.CuArray{Float64,2}:
    ...

julia> Calculators.distance_matrix_serial(pose, verlet_list)
N×N CUDA.CuArray{Float64,2}:
    ...
```
"""
function distance_matrix_serial(coords::Matrix{T}) where {T <: AbstractFloat}
    # coords must be in AoS format

    natoms = size(coords)[2]
    distance_matrix = zeros(natoms, natoms)

    for i in 1:(natoms - 1)

        for j in (i + 1):natoms
            rij = coords[:, j] - coords[:, i]
            d = sqrt(sum(rij.*rij))
            distance_matrix[i, j] = d
        end
    end
    
    return distance_matrix
end

distance_matrix_serial(pose::Pose) = distance_matrix_serial(pose.state.x)
distance_matrix_serial(state::State{T}) where {T <: AbstractFloat} = begin
    distance_matrix_serial(state.x)
end

distance_matrix_serial(pose::Pose, selection::ProtoSyn.AbstractSelection) = begin

    s = selection(pose, gather = true)
    coords = zeros(eltype(pose.state), 3, length(s))
    for (i, atom) in enumerate(s)
        coords[:, i] = pose.state[atom].t
    end

    return distance_matrix_serial(coords)
end


function distance_matrix_serial(coords::Matrix{T}, vlist::VerletList) where {T <: AbstractFloat}
    
    natoms = size(coords)[2]
    distance_matrix = zeros(natoms, natoms)
    
    @inbounds for i = 1:natoms-1

        # ptr -> location of the first neighbor of atom i
        # ptr_stop -> location of the last neighbor of atom i
        ptr = vlist.offset[i]
        if vlist.list[ptr] < 1
            continue
        end
        ptr_stop = vlist.offset[i + 1] - 2

        # load coordinates for the i-th atom
        @nexprs 3 u -> ri_u = coords[u, i]
        
        while ptr <= ptr_stop
            j = vlist.list[ptr]

            # load coordinates for the j-th atom
            # and calculate the ij vector
            @nexprs 3 u -> rij_u = coords[u, j] - ri_u
            
            # calculate the squared distance. Skip
            # if greater than cutoff
            dij_sq = @reduce 3 (+) u -> rij_u*rij_u
            distance_matrix[i, j] = sqrt(dij_sq)
            
            ptr += 1
        end
    end
    
    return distance_matrix
end

distance_matrix_serial(state::State{T}, vlist::VerletList) where {T <: AbstractFloat} = begin
    distance_matrix_serial(state.x, vlist)
end

distance_matrix_serial(pose::Pose, vlist::VerletList) = begin
    distance_matrix_serial(pose.state.x, vlist)
end