# using CUDA: sqrt, CuDeviceArray, CuDeviceMatrix, CuArray, @cuda
using CUDA

"""
    Calculators.distance_matrix_kernel(coords::CuDeviceArray{T}, distance_matrix::CuDeviceMatrix{T}, n::Int) where {T <: AbstractFloat}

Kernel for CUDA calculation of distance matrix. Should not be used directly.
Instead use `distance_matrix_cuda`.
"""
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


export distance_matrix_cuda

"""
    Calculators.distance_matrix_cuda(coords::Vector{T}) where {T <: AbstractFloat}

Return a distance_matrix with the distance of all pairs of coordinates in
`coords` (this should be a Vector{T} in AoS format), using the GPU and CUDA API.

    Calculators.distance_matrix_cuda(state::State{T})

Return a distance_matrix with the distance of all pairs of atoms in the State
`state`, using the GPU and CUDA API.

    Calculators.distance_matrix_cuda(pose::Pose)

Return a distance_matrix with the distance of all pairs of atoms in the Pose
`pose`, using the GPU and CUDA API.

    Calculators.distance_matrix_cuda(pose::Pose, selection::ProtoSyn.AbstractSelection)

Return a distance_matrix with the distance of all pairs of atoms in the Pose
`pose` who are included in the given `selection`, using the GPU and CUDA API.

# See also
`distance_matrix_serial`, `distance_matrix_simd`

# Examples
```jldoctest
julia> Calculators.distance_matrix_cuda(pose.state.x[:])
N×N CUDA.CuArray{Float64,2}:
    ...

julia> Calculators.distance_matrix_cuda(pose, an"CA")
N×N CUDA.CuArray{Float64,2}:
    ...
```
"""
function distance_matrix_cuda(coords::Vector{T}) where {T <: AbstractFloat}
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

distance_matrix_cuda(pose::Pose, selection::ProtoSyn.AbstractSelection) = begin
    
    s = selection(pose, gather = true)
    coords = zeros(eltype(pose.state), 3, length(s))
    for (i, atom) in enumerate(s)
        coords[:, i] = pose.state[atom].t
    end

    return distance_matrix_cuda(coords[:])
end

distance_matrix_cuda(pose::Pose) = distance_matrix_cuda(pose.state.x[:])
distance_matrix_cuda(state::State{T}) where {T <: AbstractFloat} = begin
    distance_matrix_cuda(state.x[:])
end