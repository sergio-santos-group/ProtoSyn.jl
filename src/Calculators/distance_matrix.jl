export distance_matrix
using CUDA
using SIMD

# ------------- SISD -----------------------------------------------------------

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

    Calculators.distance_matrix_serial(coords::Matrix{T}, verlet_list::VerletList) where {T <: AbstractFloat}
    Calculators.distance_matrix_serial(state::State{T}, verlet_list::VerletList) where {T <: AbstractFloat}
    Calculators.distance_matrix_serial(pose::Pose, verlet_list::VerletList)

Return a distance matrix with the distance of all pairs of coordinates in the
`VerletList` `verlet_list`, in the CPU using SISD architecture. Check `VerletList` for 
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
function distance_matrix(::Type{ProtoSyn.SISD_0}, coords::Matrix{T}) where {T <: AbstractFloat}
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

distance_matrix(::Type{ProtoSyn.SISD_0}, pose::Pose) = distance_matrix_serial(ProtoSyn.SISD_0, pose.state.x)
distance_matrix(::Type{ProtoSyn.SISD_0}, state::State{T}) where {T <: AbstractFloat} = begin
    distance_matrix(ProtoSyn.SISD_0, state.x)
end

distance_matrix(::Type{ProtoSyn.SISD_0}, pose::Pose, selection::ProtoSyn.AbstractSelection) = begin

    s = selection(pose, gather = true)
    coords = zeros(eltype(pose.state), 3, length(s))
    for (i, atom) in enumerate(s)
        coords[:, i] = pose.state[atom].t
    end

    return distance_matrix_serial(ProtoSyn.SISD_0, coords)
end


function distance_matrix(::Type{ProtoSyn.SISD_0}, coords::Matrix{T}, verlet_list::VerletList) where {T <: AbstractFloat}
    
    natoms = size(coords)[2]
    distance_matrix = zeros(natoms, natoms)
    
    @inbounds for i = 1:natoms-1

        # ptr -> location of the first neighbor of atom i
        # ptr_stop -> location of the last neighbor of atom i
        ptr = verlet_list.offset[i]
        if verlet_list.list[ptr] < 1
            continue
        end
        ptr_stop = verlet_list.offset[i + 1] - 2

        # load coordinates for the i-th atom
        @nexprs 3 u -> ri_u = coords[u, i]
        
        while ptr <= ptr_stop
            j = verlet_list.list[ptr]

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

distance_matrix(::Type{ProtoSyn.SISD_0}, state::State{T}, verlet_list::VerletList) where {T <: AbstractFloat} = begin
    distance_matrix(ProtoSyn.SISD_0, state.x, verlet_list)
end

distance_matrix(::Type{ProtoSyn.SISD_0}, pose::Pose, verlet_list::VerletList) = begin
    distance_matrix(ProtoSyn.SISD_0, pose.state.x, verlet_list)
end

# ------------- SIMD -----------------------------------------------------------


function distance_matrix(::Type{ProtoSyn.SIMD_1}, coords::Vector{T}) where {T <: AbstractFloat}
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

distance_matrix(::Type{ProtoSyn.SIMD_1}, pose::Pose) = distance_matrix(ProtoSyn.SIMD_1, pose.state.x[:])
distance_matrix(::Type{ProtoSyn.SIMD_1}, state::State{T}) where {T <: AbstractFloat} = begin
    distance_matrix(ProtoSyn.SIMD_1, state.x[:])
end


distance_matrix(::Type{ProtoSyn.SIMD_1}, pose::Pose, selection::ProtoSyn.AbstractSelection) = begin

    s = selection(pose, gather = true)
    coords = zeros(eltype(pose.state), 3, length(s))
    for (i, atom) in enumerate(s)
        coords[:, i] = pose.state[atom].t
    end

    return distance_matrix(ProtoSyn.SIMD_1, coords[:])
end


function distance_matrix(::Type{ProtoSyn.SIMD_1}, coords::Vector{T}, verlet_list::VerletList) where {T <: AbstractFloat}

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

distance_matrix(::Type{ProtoSyn.SIMD_1}, state::State{T}, verlet_list::VerletList) where {T <: AbstractFloat} = begin
    distance_matrix(ProtoSyn.SIMD_1, state.x[:], verlet_list)
end

distance_matrix(::Type{ProtoSyn.SIMD_1}, pose::Pose, verlet_list::VerletList) = begin
    distance_matrix(ProtoSyn.SIMD_1, pose.state.x[:], verlet_list)
end

# ------------- CUDA -----------------------------------------------------------


function distance_matrix_kernel(coords::CuDeviceArray{T}, distance_matrix::CuDeviceMatrix{T}, n::Int) where {T <: AbstractFloat}
    # Note: coords must be in AoS format

    i = ((blockIdx().y - 1) * blockDim().y) + threadIdx().y
    j = ((blockIdx().x - 1) * blockDim().x) + threadIdx().x
    
    if i <= n && j <= n

        _i = i << 2 - (2 + i)
        _j = j << 2 - (2 + j)

        dijx    = coords[_i] - coords[_j]
        dijx_sq = dijx * dijx
        dijy    = coords[_i + 1] - coords[_j + 1]
        dijy_sq = dijy * dijy
        dijz    = coords[_i + 2] - coords[_j + 2]
        dijz_sq = dijz * dijz

        dij_sq = dijx_sq + dijy_sq + dijz_sq
        distance_matrix[i, j] = CUDA.sqrt(dij_sq)
    end

    return
end


function distance_matrix(::Type{ProtoSyn.CUDA_2}, coords::Vector{T}) where {T <: AbstractFloat}
    # coords must be in AoS format
    
    _size   = trunc(Int64, length(coords)/3)
    _coords = CuArray(coords)
    
    # Define the configuration
    n_threads     = min(_size, 32)
    threads       = (n_threads, n_threads)
    n_blocks      = ceil(Int, _size / n_threads)
    blocks        = (n_blocks, n_blocks)
    
    results = CuArray(zeros(T, _size, _size))

    # println("   Size: $_size")
    # println("Threads: $threads")
    # println(" Blocks: $blocks")
    # return 
    
    @cuda blocks = blocks threads = threads distance_matrix_kernel(_coords, results, _size)
    
    return results
end

distance_matrix(::Type{ProtoSyn.CUDA_2}, pose::Pose, selection::ProtoSyn.AbstractSelection) = begin
    
    s = selection(pose, gather = true)
    coords = zeros(eltype(pose.state), 3, length(s))
    for (i, atom) in enumerate(s)
        coords[:, i] = pose.state[atom].t
    end

    return distance_matrix(ProtoSyn.CUDA_2, coords[:])
end

distance_matrix(::Type{ProtoSyn.CUDA_2}, pose::Pose) = begin
    distance_matrix(ProtoSyn.CUDA_2, pose.state.x[:])
end

distance_matrix(::Type{ProtoSyn.CUDA_2}, state::State{T}) where {T <: AbstractFloat} = begin
    distance_matrix(ProtoSyn.CUDA_2, state.x[:])
end

# ------------- DYNAMIC---------------------------------------------------------

function distance_matrix(coords::Vector{T}) where {T <: AbstractFloat}
    distance_matrix(ProtoSyn.acceleration, coords)
end

distance_matrix(pose::Pose, selection::ProtoSyn.AbstractSelection) = begin
    distance_matrix(ProtoSyn.acceleration, pose, selection)
end

distance_matrix(pose::Pose) = begin
    distance_matrix(ProtoSyn.acceleration, pose)
end

distance_matrix(state::State{T}) where {T <: AbstractFloat} = begin
    distance_matrix(ProtoSyn.acceleration, state)
end

distance_matrix(state::State{T}, verlet_list::VerletList) where {T <: AbstractFloat} = begin
    if ProtoSyn.acceleration == Type{ProtoSyn.CUDA_2}
        distance_matrix(ProtoSyn.SIMD_1, state.x, verlet_list)
    else
        distance_matrix(ProtoSyn.acceleration, state.x, verlet_list)
    end
end

distance_matrix(pose::Pose, verlet_list::VerletList) = begin
    if ProtoSyn.acceleration == Type{ProtoSyn.CUDA_2}
        distance_matrix(ProtoSyn.SIMD_1, pose, verlet_list)
    else
        distance_matrix(ProtoSyn.acceleration, pose, verlet_list)
    end
end

distance_matrix(coords::Vector{T}, verlet_list::VerletList) where {T <: AbstractFloat} = begin
    if ProtoSyn.acceleration == Type{ProtoSyn.CUDA_2}
        distance_matrix(ProtoSyn.SIMD_1, coords, verlet_list)
    else
        distance_matrix(ProtoSyn.acceleration, coords, verlet_list)
    end
end