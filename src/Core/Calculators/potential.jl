using CUDA
using LinearAlgebra

function apply_potential_kernel(coords::CuDeviceArray{T}, energies::CuDeviceArray{T}, forces::CuDeviceArray{T}, n::Int, potential::Function) where {T <: AbstractFloat}
    # ! Note: coords must be in AoS format;
    # ! Note: `potential` function receives a distance::T and returns an
    # ! energy value e::T. If it also receives a distance tuple
    # ! v::Tuple{T, T, T}, it returns the forces f1::Tuple{T, T, T} and
    # ! f2::Tuple{T, T, T} felt on both ends of the vector, based on the given
    # ! distance;

    i = ((blockIdx().y - 1) * blockDim().y) + threadIdx().y
    j = ((blockIdx().x - 1) * blockDim().x) + threadIdx().x

    if i <= n && j <= n

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
        
        # ! Note: When dealing with GPU optimizations, CUDA.jl can't allocate
        # ! new vectors. Therefore, only Tuples can be used.
        # * Note that in the next line, the returned values are the output of
        # * calling the `potential` function.
        energies[i, j], (forces[i, j, 1], forces[i, j, 2], forces[i, j, 3]), (forces[j, i, 1], forces[j, i, 2], forces[j, i, 3]) = potential(dij, v = (dijx, dijy, dijz))
    end

    return nothing
end # function

# * Coords / No Mask
"""
    apply_potential([::Type{A}], coords::Vector{T}, potential::Function) where {A <: ProtoSyn.AbstractAccelerationType, <: AbstractFloat}
    apply_potential([::Type{A}], coords::Vector{T}, potential::Function, mask::Union{ProtoSyn.Mask{C}, Matrix{T}}) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat, C <: ProtoSyn.AbstractContainer}

Apply the given `potential` to the provided `coords`, return the total energy of
the system and matrix of forces felt on each atom (forces are always calculated).
If a `mask` is provided, the resulting energy and force matrix are multiplied by
this `mask` (See [Available masks](@ref) for a list of default masks is provided
in `ProtoSyn.Calculators` module). The `potential` function should receive a
`distance::T` and return an energy value `e::T`. If it receives an optional
tuple `v::Tuple{T, T, T}`, it should also return the forces `f1::Tuple{T, T, T}`
and `f2::Tuple{T, T, T}` felt on both ends of the vector, based on the given
`distance::T` (See [Available potentials](@ref) for a list of default potential
functions available in `ProtoSyn.Calculators` module and [Creating custom
potential functions](@ref) for the correct function signatures of new
potentials). An optional parameter `Type{<: AbstractAccelerationType}` can be
provided, stating the acceleration type used to calculate this energetic
contribution (See [ProtoSyn acceleration types](@ref)).

*Selection | Mask*

    apply_potential([::Type{A}], pose::Pose, potential::Function, mask::Union{ProtoSyn.Mask{C}, Matrix{T}}, selection::AbstractSelection) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat, C <: ProtoSyn.AbstractContainer}
    apply_potential([::Type{A}], pose::Pose, potential::Function, mask::Function, selection::AbstractSelection) where {A <: ProtoSyn.AbstractAccelerationType}

*Selection | No Mask*

    apply_potential([::Type{A}], pose::Pose, potential::Function, mask::Nothing, selection::AbstractSelection) where {A <: ProtoSyn.AbstractAccelerationType}
    apply_potential([::Type{A}], pose::Pose, potential::Function, selection::AbstractSelection) where {A <: ProtoSyn.AbstractAccelerationType}

Apply the given `potential` to the selected atoms of [`Pose`](@ref) `pose` via
the provided `selection`, return the total energy of the system and matrix of
forces felt on each atom. Optionally, multiply the results by a `mask` (See
[Available masks](@ref)). If given (and not equal to `nothing`), the `mask` size
must match the `N` selected atoms. Alternatively, the given `mask` can be a
`Function`, in which case it receives a [`Pose`](@ref) `pose` as input (For the
correct signature of this `Function` `mask`, see [Creating custom masks](@ref)).

*No Selection | Mask*

    apply_potential([::Type{A}], pose::Pose, potential::Function, mask::Union{ProtoSyn.Mask{C}, Matrix{T}}) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat, C <: ProtoSyn.AbstractContainer}
    apply_potential([::Type{A}], pose::Pose, potential::Function, mask::Function) where {A <: ProtoSyn.AbstractAccelerationType}
    apply_potential([::Type{A}], pose::Pose, potential::Function, selection::Nothing, mask::Union{ProtoSyn.Mask{C}, Matrix{T}}) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat, C <: ProtoSyn.AbstractContainer}
    apply_potential([::Type{A}], pose::Pose, potential::Function, selection::Nothing, mask::Function) where {A <: ProtoSyn.AbstractAccelerationType}
    
*No Selection | No Mask*

    apply_potential([::Type{A}], pose::Pose, potential::Function, mask::Nothing) where {A <: ProtoSyn.AbstractAccelerationType}
    apply_potential([::Type{A}], pose::Pose, potential::Function) where {A <: ProtoSyn.AbstractAccelerationType}
    apply_potential([::Type{A}], pose::Pose, potential::Function, selection::Nothing, mask::Nothing) where {A <: ProtoSyn.AbstractAccelerationType}
    apply_potential([::Type{A}], pose::Pose, potential::Function, selection::Nothing) where {A <: ProtoSyn.AbstractAccelerationType}

Apply the given `potential` to the all atoms of [`Pose`](@ref) `pose`, return
the total energy of the system and matrix of forces felt on each [`Atom`](@ref).
If given (and not equal to `nothing`), the `mask` size must match the total
number of [`Atom`](@ref) instances in the pose. Alternatively, the given `mask`
can be a`Function`, in which case it receives a [`Pose`](@ref) `pose` as input
(For the correct signature of this `Function` `mask`, see
[Creating custom masks](@ref)).

!!! ukw "Note:"
    As of ProtoSyn 1.0, this function's acceleration type defaults to
    `CUDA_2` regardless of the requested acceleration type. This may be changed
    in future iterations.

# Examples
```jldoctest
julia> fbr = ProtoSyn.Calculators.get_flat_bottom_potential(2.0, 5.0)
...

julia> sidechain = !an"^CA\$|^N\$|^C\$|^H\$|^O\$"r
UnarySelection{ProtoSyn.Stateless}(!, FieldSelection{ProtoSyn.Stateless,Atom}(r"^CA\$|^N\$|^C\$|^H\$|^O\$", :name, occursin))

julia> mask = ProtoSyn.Calculators.intra_residue_mask(pose, sidechain)

julia> e, f = ProtoSyn.Calculators.apply_potential(ProtoSyn.CUDA_2, pose, fbr, mask, sidechain)
(2.1792609341377363, [11.380 … -74.232])
```
"""
function apply_potential(::Type{ProtoSyn.CUDA_2}, coords::Vector{T}, potential::Function) where {T <: AbstractFloat}

    _size   = trunc(Int64, length(coords)/3)
    _coords = CuArray(coords)
    
    # Define the configuration
    n_threads     = min(_size, 20) # ?
    threads       = (n_threads, n_threads)
    n_blocks      = ceil(Int, _size / n_threads)
    blocks        = (n_blocks, n_blocks)

    forces   = CuArray(zeros(T, _size, _size, 3))
    energies = CuArray(zeros(T, _size, _size))

    @cuda blocks = blocks threads = threads apply_potential_kernel(_coords, energies, forces, _size, potential)

    return energies, forces
end # function

apply_potential(::Type{A}, coords::Vector{T}, potential::Function) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat} = begin
    apply_potential(ProtoSyn.CUDA_2, coords, potential)
end


# * Coords / Mask
function apply_potential(::Type{ProtoSyn.CUDA_2}, coords::Vector{T}, potential::Function, mask::Union{ProtoSyn.Mask{C}, Matrix{T}}) where {T <: AbstractFloat, C <: ProtoSyn.AbstractContainer}

    e, f   = apply_potential(ProtoSyn.CUDA_2, coords, potential)

    @assert size(e) == size(mask) "Tried to apply a mask but the size $(size(mask.content)) doesn't match with the input size $(size(e))."

    energy = sum(e .* CuArray(mask))
    
    f = map(*, f, CuArray(repeat(mask, outer = (1, 1, 3))))
    forces = collect(reshape(sum(f, dims = 1), size(f)[1], 3)')

    return energy, forces
end # function

apply_potential(::Type{A}, coords::Vector{T}, potential::Function, mask::Union{ProtoSyn.Mask{C}, Matrix{T}}) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat, C <: ProtoSyn.AbstractContainer} = begin
    apply_potential(ProtoSyn.CUDA_2, coords, potential, mask)
end

# ---


# * Selection / No Mask
apply_potential(::Type{ProtoSyn.CUDA_2}, pose::Pose, potential::Function, selection::AbstractSelection) = begin
    sele = ProtoSyn.promote(selection, Atom)(pose).content
    coords = pose.state.x.coords[:, sele][:]

    e, _f = apply_potential(ProtoSyn.CUDA_2, coords, potential)
    energy = sum(e)
    _f = collect(reshape(sum(_f, dims = 1), size(_f)[1], 3)')

    forces = zeros(eltype(pose.state), 3, pose.state.size)
    forces[:, sele] = _f
    return energy, forces
end

apply_potential(::Type{ProtoSyn.CUDA_2}, pose::Pose, potential::Function, mask::Nothing, selection::AbstractSelection) = begin
    apply_potential(ProtoSyn.CUDA_2, pose, potential, selection)
end

apply_potential(::Type{A}, pose::Pose, potential::Function, selection::AbstractSelection) where {A <: ProtoSyn.AbstractAccelerationType} = begin
    apply_potential(ProtoSyn.CUDA_2, pose, potential, selection)
end

apply_potential(::Type{A}, pose::Pose, potential::Function, mask::Nothing, selection::AbstractSelection) where {A <: ProtoSyn.AbstractAccelerationType} = begin
    apply_potential(ProtoSyn.CUDA_2, pose, potential, selection)
end

# ---

# * Selection / Mask
apply_potential(::Type{ProtoSyn.CUDA_2}, pose::Pose, potential::Function, mask::Union{ProtoSyn.Mask{C}, Matrix{T}}, selection::AbstractSelection) where {T <: AbstractFloat, C <: ProtoSyn.AbstractContainer} = begin
    sele = ProtoSyn.promote(selection, Atom)(pose).content
    coords = pose.state.x.coords[:, sele][:]
    e, _f = apply_potential(ProtoSyn.CUDA_2, coords, potential, mask)
    f = zeros(eltype(pose.state), 3, pose.state.size)
    f[:, sele] = _f
    return e, f
end

apply_potential(::Type{ProtoSyn.CUDA_2}, pose::Pose, potential::Function, mask::Function, selection::AbstractSelection) = begin
    apply_potential(ProtoSyn.CUDA_2, pose, potential, mask(pose), selection)
end

apply_potential(::Type{A}, pose::Pose, potential::Function, mask::Union{ProtoSyn.Mask{C}, Matrix{T}}, selection::AbstractSelection) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat, C <: ProtoSyn.AbstractContainer} = begin
    apply_potential(ProtoSyn.CUDA_2, pose, potential, mask, selection)
end

apply_potential(::Type{A}, pose::Pose, potential::Function, mask::Function, selection::AbstractSelection) where {A <: ProtoSyn.AbstractAccelerationType} = begin
    apply_potential(ProtoSyn.CUDA_2, pose, potential, mask, selection)
end

# ---


# * No Selection / Mask
apply_potential(::Type{A}, pose::Pose, potential::Function, mask::Union{ProtoSyn.Mask{C}, Matrix{T}}) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat, C <: ProtoSyn.AbstractContainer} = begin
    coords = pose.state.x.coords[:]
    apply_potential(ProtoSyn.CUDA_2, coords, potential, mask)
end

apply_potential(::Type{A}, pose::Pose, potential::Function, selection::Nothing, mask::Union{ProtoSyn.Mask{C}, Matrix{T}}) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat, C <: ProtoSyn.AbstractContainer} = begin
    coords = pose.state.x.coords[:]
    apply_potential(ProtoSyn.CUDA_2, coords, potential, mask)
end

apply_potential(::Type{A}, pose::Pose, potential::Function, mask::Function) where {A <: ProtoSyn.AbstractAccelerationType} = begin
    apply_potential(ProtoSyn.CUDA_2, pose, potential, mask(pose))
end

apply_potential(::Type{A}, pose::Pose, potential::Function, selection::Nothing, mask::Function) where {A <: ProtoSyn.AbstractAccelerationType} = begin
    apply_potential(ProtoSyn.CUDA_2, pose, potential, mask(pose))
end


# * No Selection / No Mask
apply_potential(::Type{A}, pose::Pose, potential::Function) where {A <: ProtoSyn.AbstractAccelerationType} = begin
    coords = pose.state.x.coords[:]
    e, _f = apply_potential(ProtoSyn.CUDA_2, coords, potential)
    energy = sum(e)
    forces = collect(reshape(sum(_f, dims = 1), size(_f)[1], 3)')
    return energy, forces
end

apply_potential(::Type{A}, pose::Pose, potential::Function, selection::Nothing) where {A <: ProtoSyn.AbstractAccelerationType} = begin
    apply_potential(ProtoSyn.CUDA_2, pose, potential)
end

apply_potential(::Type{A}, pose::Pose, potential::Function, selection::Nothing, mask::Nothing) where {A <: ProtoSyn.AbstractAccelerationType} = begin
    apply_potential(ProtoSyn.CUDA_2, pose, potential)
end



# ------------------------------------------------------------------------------
# * Available potential functions

"""
    get_flat_bottom_potential(;d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf) where {T <: AbstractFloat}

Return a flat-bottom potential function, using the specified distances. The
potential is made up of 5 different sectors, each with the following functions:

\$f_{1}) \\;\\;\\;\\;\\;\\; e = m_{1} \\cdot d + b_{1} \\;\\;\\;\\;\\;\\; \\left \\{ d < d_{1} \\right \\}\\,\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

\$f_{2}) \\;\\;\\;\\;\\;\\; e = \\left (d-d_{2}  \\right )^{2} \\;\\;\\;\\;\\;\\;\\;\\; \\left \\{ d_{1} \\leqslant d < d_{2} \\right \\}\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

\$f_{3}) \\;\\;\\;\\;\\;\\; e = 0 \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\; \\left \\{ d_{2} \\leqslant d \\leqslant d_{3} \\right \\}\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

\$f_{4}) \\;\\;\\;\\;\\;\\; e = \\left (d-d_{3}  \\right )^{2} \\;\\;\\;\\;\\;\\;\\;\\; \\left \\{ d_{4} < d \\leqslant d_{4} \\right \\}\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

\$f_{5}) \\;\\;\\;\\;\\;\\; e = m_{2} \\cdot d + b_{2} \\;\\;\\;\\;\\;\\;  \\left \\{ d > d_{4} \\right \\}\\,\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

Where 

\$m_{1} = 2 \\left ( d_{1}-d_{2} \\right ) \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$
\$b_{1} = f_{2}\\left ( d_{1} \\right ) - m_{1} \\cdot d_{1} \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$
\$m_{2} = 2\\left ( d_{4} - d_{3} \\right ) \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$
\$b_{2} = f_{4}\\left ( d_{4} \\right ) - m_{2} \\cdot d_{4} \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

*The resulting function can be called with the following signature:*
    
```jldoctest
flat_bottom_potential(d::T; v::Opt{Vector{T}} = nothing) where {T <: AbstractFloat}
```

Return an energy value based on the provided distance `d`. If a vector `v` is
also provided (optional), the flat-bottom restraint will also return the forces
`f1` and `f2` (the forces felt on both ends of the vector `v`). The vector `v`
should have length = 3, corresponding to the 3 dimensions of the distance
between the two [`Atom`](@ref) instances (X, Y and Z). For more information on
the flat-bottom potential, please read:
[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4692055/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4692055/). 

# See also
[`apply_potential`](@ref ProtoSyn.Calculators.apply_potential) [`calc_flat_bottom_restraint`](@ref ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint)

# Examples
```jldoctest
julia> f = Calculators.get_flat_bottom_potential(1.0, 2.0, 3.0, 4.0)

julia> f(2.5)
0.0

julia> f(1.73, v = [1.0, 1.0, 1.0])
(0.0729, [-0.54, -0.54, -0.54], [0.54, 0.54, 0.54])
```
"""
function get_flat_bottom_potential(;d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf) where {T <: AbstractFloat}
    
    @assert d2 >= d1 "d2 must be of equal value or over d1"
    @assert d3 >= d2 "d3 must be of equal value or over d2"
    @assert d4 >= d3 "d4 must be of equal value or over d3"

    return function flat_bottom_potential(d::T; v::Opt{Tuple{T, T, T}} = nothing) where {T <: AbstractFloat}

        v !== nothing && begin
            f1 = (T(0.0), T(0.0), T(0.0))
            f2 = (T(0.0), T(0.0), T(0.0))
        end

        
        d12     = d1 - d2
        m1      = 2 * d12
        b1      = d12 * d12 - m1 * d1
        
        d43     = d4 - d3
        m2      = 2 * d43
        b2      = d43 * d43 - m2 * d4
        
        if d < d1 # * Linear Left (1)
            e = m1 * d + b1
            
            v !== nothing && begin
                f1 = .- v .* m1
                f2 =    v .* m1
            end
        elseif d1 <= d < d2 # * Quadratic Left (2)
            e = (d - d2) * (d - d2)
            
            v !== nothing && begin
                    factor1 = 2 * (d - d2)
                    f1      = .- v .* factor1
                    f2      =    v .* factor1
                end
            elseif d2 <= d <= d3 # * Flat (3)
                e = 0
            elseif d3 < d <= d4 # * Quadratic Right (4)
                e = (d - d3) * (d - d3)
                
                v !== nothing && begin
                factor2 = 2 * (d - d3)
                f1      = .- v .* factor2
                f2      =    v .* factor2
            end            
        else # * Linear Right (5)
            e = m2 * d + b2
            
            v !== nothing && begin
            f1 = .- v .* m2
            f2 =    v .* m2
        end
    end
    
    v !== nothing && return e, f1, f2
    return e
    end
end

# ------------------------------------------------------------------------------
# * Available potential masks

"""
    intra_residue_mask(pose::Pose, selection::AbstractSelection)

For all the atoms in the provided `AbstractSelection` `selection` (N), return a
2D N x N [`Mask`](@ref) with all the atoms of the given [`Pose`](@ref) `pose`
not in the same residue selected. 

!!! ukw "Note:"
    This function is rather heavy and has low performance. If no design effort
    is being made (where the sequence changes), the resulting [`Mask`](@ref)
    from this function can and should be re-used (only calculated once). If, for
    a specific application, the `AbstractSelection` `selection` remains constant
    but the [`Mask`](@ref) needs to be re-calculated (for example, because there
    was a design/mutation step, use the _functor_ resulting from
    [`get_intra_residue_mask`](@ref)).

# See also
[`diagonal_mask`](@ref) [`get_intra_residue_mask`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Calculators.intra_residue_mask(pose, !an"^CA\$|^N\$|^C\$|^H\$|^O\$"r)
ProtoSyn.Mask{Atom}(Bool[1 1 … 0 0; 1 1 … 0 0; … ; 0 0 … 1 1; 0 0 … 1 1])
```
"""
function intra_residue_mask(pose::Pose, selection::AbstractSelection)
    # ! Note: This function is rather heavy. If no design effort is being made
    # ! (where the sequence changes), the resulting map from this function can
    # ! and should be re-used (only calculated once).

    s = selection(pose)
    N = count(s)
    mask = ProtoSyn.Mask{Atom}(N, N)
    for r in eachresidue(pose.graph)
        rm1 = ProtoSyn.promote(SerialSelection{Residue}(r.id, :id), Atom)(pose)
        rm1.content = rm1.content[s.content]
        mask |= ProtoSyn.cross2d(rm1)
    end

    return !mask
end


"""
    get_intra_residue_mask(selection::AbstractSelection)

Provides the [`intra_residue_mask`](@ref) function as a _functor_, which will
calculate the intra residue mask for the given `AbstractSelection` `selection`.
Useful when creating a new [`EnergyFunctionComponent`](@ref) or when the
[`Mask`](@ref) should be updated each step/call.

# See also
[`intra_residue_mask`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Calculators.get_intra_residue_mask(!an"^CA\$|^N\$|^C\$|^H\$|^O\$"r)
(::ProtoSyn.Calculators.var"#_intra_residue_mask#5"{UnarySelection{ProtoSyn.Stateless}}) (generic function with 1 method)
```
"""
function get_intra_residue_mask(selection::AbstractSelection)
    return function _intra_residue_mask(pose::Pose)
        s = selection(pose)
        N = count(s)
        mask = ProtoSyn.Mask{Atom}(N, N)
        for r in eachresidue(pose.graph)
            rm1 = ProtoSyn.promote(SerialSelection{Residue}(r.id, :id), Atom)(pose)
            rm1.content = rm1.content[s.content]
            mask |= ProtoSyn.cross2d(rm1)
        end

        return !mask
    end
end


"""
    diagonal_mask(pose::Pose, selection::AbstractSelection)

For all the atoms in the provided `AbstractSelection` `selection` (N), return a
2D N x N [`Mask`](@ref) with all the [`Atom`](@ref) instances of the given
[`Pose`](@ref) `pose` not in the natural diagonal selected (i.e. ignores same
atom interaction artifacts).

!!! ukw "Note:"
    When the selection is constant but the resulting [`Mask`](@ref) needs to be
    re-calculated every call/step (for example, due to a design or mutation
    step), consider using the _functor_ from [`get_diagonal_mask`](@ref).

# See also
[`intra_residue_mask`](@ref) [`get_diagonal_mask`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Calculators.diagonal_mask(pose, an"CA")
ProtoSyn.Mask{Atom}(3, 3)
3×3 BitArray{2}:
 0  1  1
 1  0  1
 1  1  0
```
"""
function diagonal_mask(pose::Pose, selection::AbstractSelection)

    N = count(selection(pose))
    return !ProtoSyn.Mask{Atom}(BitArray(Matrix{Bool}(I, N, N)))
end


"""
    get_diagonal_mask(selection::AbstractSelection)

Provides the [`diagonal_mask`](@ref) as a _functor_, which will calculate the
diagonal mask for the given `AbstractSelection` `selection`. Useful when
creating a new [`EnergyFunctionComponent`](@ref) or when the [`Mask`](@ref)
should be updated each step/call.

# See also
[`diagonal_mask`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Calculators.get_diagonal_mask(an"CA")
(::ProtoSyn.Calculators.var"#_diagonal_mask#6"{FieldSelection{ProtoSyn.Stateless,Atom}}) (generic function with 1 method)
```
"""
function get_diagonal_mask(selection::AbstractSelection)
    return function _diagonal_mask(pose::Pose)
        N = count(selection(pose))
        return !ProtoSyn.Mask{Atom}(BitArray(Matrix{Bool}(I, N, N)))
    end
end


"""
    load_map([::Type{T}], filename::String) where {T <: AbstractFloat}

Load the map in the `filename` file (i.e. Contact Map). The file should be in
PFRMAT RR format (See: [https://predictioncenter.org/casp13/index.cgi?page=format#RR](https://predictioncenter.org/casp13/index.cgi?page=format#RR)).
Returns an N x N map of the found weights, with pairs not identified in the file
set to 0.0 (N is the maximum indentifier found on the file. As an example, it
might be the case where a peptide has 74 residues, but no pair with residue 74
is found on the file, the maximum identifier found might be 72, for example. In
this case, the resulting map will have size 72 x 72. In order to ensure the
loaded map size matches the underlying peptide size, consider adding an entry of
0.0 on the map file, with the correct maximum identifier). Note: If no
optional type `T` is provided, will use `ProtoSyn.Units.defaultFloat`.

# Examples
```jldoctest
julia> cmap = ProtoSyn.Calculators.load_map("contact_map_example.txt")
73×73 Array{Float64,2}:
 ...
```
"""
function load_map(::Type{T}, filename::String) where {T <: AbstractFloat}

    _map = Dict{Tuple{Int, Int}, T}()

    open(filename, "r") do map_file
        for line in eachline(map_file)
            elems = split(line)
            length(elems[1]) > 4 && continue
            elems[1] == "END" && continue
            try
                α = parse(T, elems[5])
                _map[(parse(Int, elems[1]), parse(Int, elems[2]))] = α
            catch BoundsError
                continue
            end
        end
    end

    k    = keys(_map)
    N    = maximum((maximum(last, k), maximum(first, k)))
    map  = zeros((N, N))

    for (key, value) in _map
        map[key[1], key[2]] = value
        map[key[2], key[1]] = value
    end

    return map
end

load_map(filename::String) = load_map(ProtoSyn.Units.defaultFloat, filename)