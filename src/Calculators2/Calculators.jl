module Calculators

using ..ProtoSyn
using Base.Cartesian

include("verlet_list.jl")
include("serial.jl")

@info "Loading SIMD"
using SIMD
include("simd.jl")

@info "Loading CUDA"
using CUDA
include("cuda.jl")

# function distance_matrix_kernel(coords::CuDeviceVector{T}, distance_matrix::CuDeviceMatrix{T}, n::Int) where {T <: AbstractFloat}
#     # Note: Example coords are is SoA format

#     thread_id          = threadIdx().x
#     n_threads_in_block = blockDim().x
#     i                  = blockIdx().y
#     j                  = (blockIdx().x - 1) * n_threads_in_block + thread_id

#     cache = @cuDynamicSharedMem(T, 512)
    
#     if i <= n && j <= n && i < j
#         # 1. Calculate distance
#         dijx    = coords[i] - coords[j]
#         dijx_sq = dijx * dijx
#         dijy    = coords[i+n] - coords[j+n]
#         dijy_sq = dijy * dijy
#         dijz    = coords[i+2*n] - coords[j+2*n]
#         dijz_sq = dijz * dijz
#         dij_sq  = dijx_sq + dijy_sq + dijz_sq
#         cache[thread_id] = srqt(dij_sq)
#     else
#         cache[thread_id] = 0.0
#     end

#     # 4. Calculate sum of energies (in cache) for this block
#     sync_threads()
    
#     level::Int = 512 ÷ 2
#     # Note 1 ÷ 2 is the integer division, and is different from 1 / 2.

#     while level != 0
#         if thread_id <= level
#             cache[thread_id] += cache[thread_id + level]
#         end
#         sync_threads()
#         level = level ÷ 2
#     end

#     # 5. Export energy sum result
#     distance_matrix[blockIdx().x, blockIdx().y] = cache[1]

#     return
# end


# function cuda(state::State{F, T}) where {F <: AbstractMatrix, T <: AbstractFloat}
    
#     coords = CuArray(view(state.coords.values, :)) # Test is in SoA
    
#     # Define the configuration
#     total_threads = min(state.n, 512)
#     threads       = (total_threads, )
#     blocks        = (ceil(Int, state.n / total_threads), state.n)
#     shmem         = 512 * sizeof(T)
    
#     distance_matrix = CuArray(zeros(T, blocks))
    
#     @cuda blocks = blocks threads = threads shmem = shmem cuda_kernel(coords, distance_matrix, state.size)
#     return distance_matrix
# end


end