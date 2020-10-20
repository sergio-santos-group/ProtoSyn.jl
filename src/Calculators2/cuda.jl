function distance_matrix_kernel(coords::CuDeviceArray{T}, distance_matrix::CuDeviceMatrix{T}, n::Int) where {T <: AbstractFloat}
    # Note: coords must be in AoS format

    i = ((blockIdx().y - 1) * blockDim().y) + threadIdx().y
    j = ((blockIdx().x - 1) * blockDim().x) + threadIdx().x
    
    if i <= n && j <= n && i < j

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


function cuda(coords::Vector{T}) where {T <: AbstractFloat}
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

cuda(pose::Pose, selection::Opt{ProtoSyn.AbstractSelection}) = begin
    
    s = selection(pose, gather = true)
    coords = zeros(eltype(pose.state), 3, length(s))
    for (i, atom) in enumerate(s)
        coords[:, i] = pose.state[atom].t
    end

    return cuda(coords[:])
end

cuda(pose::Pose) = cuda(pose.state.x[:])
cuda(state::State{T}) where {T <: AbstractFloat} = cuda(state.x[:])