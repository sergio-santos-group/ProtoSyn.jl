# * CUDA_2

# Kernel function
function resolve_calculation_kernel(coords::CuDeviceArray{T}, energies::CuDeviceArray{T}, forces::CuDeviceArray{T}, n::Int, potential::Function, qs::CuDeviceArray{T}) where {T <: AbstractFloat}
    # ! Note: coords must be in AoS format;
    # ! Note: `potential` function receives a distance::T and returns an
    # ! energy value e::T. If it also receives a distance tuple
    # ! v::Tuple{T, T, T}, it returns the forces f1::Tuple{T, T, T} and
    # ! f2::Tuple{T, T, T} felt on both ends of the vector, based on the given
    # ! distance;

    i = ((blockIdx().y - 1) * blockDim().y) + threadIdx().y
    j = ((blockIdx().x - 1) * blockDim().x) + threadIdx().x

    if i <= n && j <= n

        if j <= i
            energies[j, i] = T(0.0)
            return nothing
        end
        _i = i << 2 - (2 + i)
        _j = j << 2 - (2 + j)

        dijx    = coords[_j] - coords[_i]
        dijx_sq = dijx * dijx
        dijy    = coords[_j + 1] - coords[_i + 1]
        dijy_sq = dijy * dijy
        dijz    = coords[_j + 2] - coords[_i + 2]
        dijz_sq = dijz * dijz

        dij_sq  = dijx_sq + dijy_sq + dijz_sq
        dij     = CUDA.sqrt(dij_sq)
        v = (dijx/dij, dijy/dij, dijz/dij) # normalization
        
        # ! Note: When dealing with GPU optimizations, CUDA.jl can't allocate
        # ! new vectors. Therefore, only Tuples can be used.
        # * Note that in the next line, the returned values are the output of
        # * calling the `potential` function.
        energies[j, i], (forces[j, i, 1], forces[j, i, 2], forces[j, i, 3]), (forces[i, j, 1], forces[i, j, 2], forces[i, j, 3]) = potential(dij, v = v, qi = qs[i], qj = qs[j])
    end

    return nothing
end

# Main function
function resolve_calculation(::Type{ProtoSyn.CUDA_2},
    pose::Pose,
    potential::Function,
    update_forces::Bool,
    verlet_list::Union{VerletList, Nothing},
    coords::Vector{T},
    mask::MaskMap,
    indexes::Vector{Int}) where {T <: AbstractFloat}

    # Note: verlet_list and update_forces are ignored in CUDA_2 mode
    # (in version 1.01)

    _size   = trunc(Int64, length(coords)/3)
    _coords = CuArray(coords)
    
    # Define the configuration
    n_threads     = min(_size, 20) # ?
    threads       = (n_threads, n_threads)
    n_blocks      = ceil(Int, _size / n_threads)
    blocks        = (n_blocks, n_blocks)

    forces   = CuArray(zeros(T, _size, _size, 3))
    energies = CuArray(zeros(T, _size, _size))

    qs = CuArray([pose.state[index].δ for index in indexes])

    @cuda blocks = blocks threads = threads resolve_calculation_kernel(_coords, energies, forces, _size, potential, qs)

    return resolve_mask(ProtoSyn.CUDA_2, pose, energies, forces, update_forces, mask)
end