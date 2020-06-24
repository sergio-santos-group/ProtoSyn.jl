using CUDA


function cuda_kernel(coords::CuDeviceVector{T}, energy::CuDeviceMatrix{T}, n::Int) where {T <: AbstractFloat}
    # Note: Example coords are is SoA format

    thread_id          = threadIdx().x
    n_threads_in_block = blockDim().x
    i                  = blockIdx().y
    j                  = (blockIdx().x - 1) * n_threads_in_block + thread_id
    σ                  = T(1)
    ϵ                  = T(1)

    cache = @cuDynamicSharedMem(T, 512)
    
    if i <= n && j <= n && i < j
        # 1. Calculate distance
        dijx    = coords[i] - coords[j]
        dijx_sq = dijx * dijx
        dijy    = coords[i+n] - coords[j+n]
        dijy_sq = dijy * dijy
        dijz    = coords[i+2*n] - coords[j+2*n]
        dijz_sq = dijz * dijz
        dij_sq  = dijx_sq + dijy_sq + dijz_sq
        
        # 2. Calculate energy
        lj2 = σ / dij_sq
        lj6 = lj2*lj2*lj2
        e   = ϵ*(lj6*lj6 - lj6)

        # 3. Save this result to the cache
        cache[thread_id] = e
    else
        cache[thread_id] = 0.0
    end

    # 4. Calculate sum of energies (in cache) for this block
    sync_threads()
    
    level::Int = 512 ÷ 2
    # Note 1 ÷ 2 is the integer division, and is different from 1 / 2.

    while level != 0
        if thread_id <= level
            cache[thread_id] += cache[thread_id + level]
        end
        sync_threads()
        level = level ÷ 2
    end

    # 5. Export energy sum result
    energy[blockIdx().x, blockIdx().y] = cache[1]

    return
end


function cuda(state::State{F, T}) where {F <: AbstractMatrix, T <: AbstractFloat}
    
    coords = CuArray(view(state.coords.values, :)) # Test is in SoA
    
    # Define the configuration
    dev           = device()
    total_threads = min(state.n, 512)
    threads       = (total_threads, )
    blocks        = (ceil(Int, state.n / total_threads), state.n)
    shmem         = 512 * sizeof(T)
    
    distances = CuArray(zeros(T, blocks))
    
    @cuda blocks = blocks threads = threads shmem = shmem cuda_kernel(coords, distances, state.n)
    return Matrix(distances)
end