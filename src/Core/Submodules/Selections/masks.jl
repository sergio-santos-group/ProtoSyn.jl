export Mask

"""
    ProtoSyn.Mask{T}(content::BitVector) where {T <: AbstractContainer}

A `Mask` is a typed wrapped on a `BitArray`, and therefore can be N dimensional.
Utility functions are provided for [`Mask`](@ref) instances up to 2D.

    ProtoSyn.Mask{T}() where {T <: AbstractContainer}

Return an empty [`Mask`](@ref) of length 0.

    ProtoSyn.Mask{T}(n::Int) where {T <: AbstractContainer}

Return a [`Mask`](@ref) of length `n` with all bits set to `false`.

    ProtoSyn.Mask{T}(n::Tuple{Int, Int}) where {T <: AbstractContainer}
    ProtoSyn.Mask{T}(n1::Int, n2::Int) where {T <: AbstractContainer}

Return a 2D [`Mask`](@ref) of size `n` (or size `(n1, n2)`) with all bits set to
`false`.

# Selection type

The selection type of a [`Mask`](@ref) instance can be any
`T <: AbstractContainer`.

!!! ukw "Note:"
    Combinatory operators (`&` and `|`) and unary operators (`!`) can be
    employed with [`Mask`](@ref) instances.

# Examples
```
julia> ProtoSyn.Mask{Atom}()
ProtoSyn.Mask{Atom}(0,)
0-element BitArray{1}

julia> ProtoSyn.Mask{Residue}((3, 3))
ProtoSyn.Mask{Residue}(3, 3)
3×3 BitMatrix:
 0  0  0
 0  0  0
 0  0  0

julia> !ProtoSyn.Mask{Atom}(5)
ProtoSyn.Mask{Atom}(5,)
5-element BitVector:
 1
 1
 1
 1
 1
```
"""
mutable struct Mask{T <: AbstractContainer}
    content::BitArray
end

Mask{T}() where {T <: AbstractContainer} = Mask{T}(BitVector())
Mask{T}(n::Int) where {T <: AbstractContainer} = Mask{T}(falses(n))
Mask{T}(n::Tuple{Int, Int}) where {T <: AbstractContainer} = Mask{T}(falses(n))
Mask{T}(n1::Int, n2::Int) where {T <: AbstractContainer} = begin
    Mask{T}(falses((n1, n2)))
end

selection_type(::Mask{T}) where {T} = T


# --- Utilities and Overloads --------------------------------------------------
Base.getindex(m::Mask{T}, i::Int) where {T <: AbstractContainer} = m.content[i]
Base.getindex(m::Mask{T}, i::UnitRange{<: Real}) where {T <: AbstractContainer} = m.content[i]

Base.setindex!(m::Mask{T}, v::Bool, i::Real) where {T <: AbstractContainer} = m.content[i] = v
Base.setindex!(m::Mask{T}, v::Bool, i::UnitRange{<: Real}) where {T <: AbstractContainer} = m.content[i] .= v

@inline Base.length(m::Mask{T}) where {T <: AbstractContainer} = length(m.content)

Base.iterate(m::Mask{T}, (s,)=(1,)) where {T <: AbstractContainer} = begin
    if s > length(m)
        return nothing
    end
    (m.content[s], (s+1,))
end

Base.copy(m::Mask{T}) where {T <: AbstractContainer} = begin
    return Mask{T}(copy(m.content))
end

Base.size(m::Mask{T}) where {T <: AbstractContainer} = begin
    return size(m.content)
end

function Base.show(io::IO, m::Mask{T}) where {T <: AbstractContainer}
    c = Base.repr(m.content, context = IOContext(stdout, :limit => true))[5:end]
    print(io, "ProtoSyn.Mask\n")
    print(io, " ├── Type: $T\n")
    print(io, " ├── Size: $(size(m.content))\n")
    print(io, " ├── Count: $(count(m.content)) / $(length(m.content))\n")
    print(io, " └── Content: $c")
end

# ---

Base.repeat(m::Mask{T}; outer::Opt{Tuple{Int, Int, Int}} = nothing) where {T <: AbstractContainer} = begin
    repeat(m.content, outer = outer)
end

CUDA.CuArray(m::Mask{T}) where {T <: AbstractContainer} = begin
    CUDA.CuArray(m.content)
end

# --- Unary operations ---------------------------------------------------------
Base.:|(m1::Mask{T}, m2::Mask{T}) where {T <: AbstractContainer} = begin
    return Mask{T}(m1.content .| m2.content)
end

Base.:&(m1::Mask{T}, m2::Mask{T}) where {T <: AbstractContainer} = begin
    return Mask{T}(m1.content .& m2.content)
end

Base.:!(m1::Mask{T}) where {T <: AbstractContainer} = begin
    m1.content = .!m1.content
    return m1
end

# --- 2D Masks utilities -------------------------------------------------------

"""
    tile!(m::Mask{T}) where {T <: AbstractContainer}

Repeat the content of the provided [`Mask`](@ref) `m` N times, where N is the
length of `m`. Note that this function expects a 1D [`Mask`](@ref) as input,
otherwise will just return the received [`Mask`](@ref). Overwrite the contents
of the input [`Mask`](@ref) with the output.

# See also
[`tile`](@ref)

# Examples
```
julia> m = ProtoSyn.Mask{Atom}(4); m[1] = true; m
ProtoSyn.Mask{Atom}(4,)
4-element BitArray{1}:
 1
 0
 0
 0

julia> ProtoSyn.tile!(m)
ProtoSyn.Mask{Atom}(4, 4)
4×4 BitArray{2}:
 1  1  1  1
 0  0  0  0
 0  0  0  0
 0  0  0  0
```
"""
function tile!(m::Mask{T}) where {T <: AbstractContainer}
    _size = size(m.content)
    if length(_size) != 1 # Case 2D matrix is received
        return
    end
    m.content = hcat(collect(Iterators.repeated(m.content, _size[1]))...)
end


"""
    tile(m::Mask{T}) where {T <: AbstractContainer}

Repeat and return the content of the provided [`Mask`](@ref) `m` N times, where
N is the length of `m`, along the horizontal dimension. Note that this function
expects a 1D [`Mask`](@ref) as input, otherwise will just return the received
[`Mask`](@ref).

# See also
[`tile!`](@ref)

# Examples
```
julia> m = ProtoSyn.Mask{Atom}(4); m[1] = true; m
ProtoSyn.Mask{Atom}(4,)
4-element BitArray{1}:
 1
 0
 0
 0

julia> ProtoSyn.tile(m)
ProtoSyn.Mask{Atom}(4, 4)
4×4 BitArray{2}:
 1  1  1  1
 0  0  0  0
 0  0  0  0
 0  0  0  0
```
"""
function tile(m::Mask{T}) where {T <: AbstractContainer}
    _size = size(m.content)
    if length(_size) != 1 # Case 2D matrix is received
        return m
    end
    return Mask{T}(hcat(collect(Iterators.repeated(m.content, _size[1]))...))
end


Base.adjoint(m::Mask{T}) where {T <: AbstractContainer} = begin
    # ! Note: This is a copy of the original mask.
    # ? Should be view?
    return Mask{T}(m.content')
end


"""
    cross2d(m1::Mask{T}, m2::Mask{T}) where {T <: AbstractContainer}

Combine both the provided [`Mask`](@ref) instances `m1` and `m2`, so that the
returned result is a 2D Mask where each position (N, M) is the result of
combining `m1[N] & m2[M]`.

    cross2d(m1::Mask{T}) where {T <: AbstractContainer}

If only 1 [`Mask`](@ref) `m1` is provided, return the `cross2d` result of `m1`
and `m1'` (`m1` transpose).

# See also
[`tile`](@ref)

# Examples
```
julia> m = ProtoSyn.Mask{Atom}(4); m[1] = m[3] = true; m
ProtoSyn.Mask{Atom}(4,)
4-element BitArray{1}:
 1
 0
 1
 0

julia> ProtoSyn.cross2d(m)
ProtoSyn.Mask{Atom}(4, 4)
4×4 BitArray{2}:
 1  0  1  0
 0  0  0  0
 1  0  1  0
 0  0  0  0
```
"""
function cross2d(m1::Mask{T}, m2::Mask{T}) where {T <: AbstractContainer}
    _m1 = tile(m1)
    _m2 = tile(m2)
    return _m1 & _m2
end

cross2d(m1::Mask{T}) where {T <: AbstractContainer} = cross2d(m1, m1')

# # --- Dynamic masks

# mutable struct DynamicMask
#     f::Function
#     sele::AbstractSelection
# end

# (dm::dynamicMask)(pose::Pose) = begin
#     return dm.f(pose, dm.sele)
# end