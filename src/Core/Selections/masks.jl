"""
    ProtoSyn.Mask{T}(content::BitVector) where {T <: AbstractContainer}

A `Mask` is a typed wrapped on a `BitArray`. `Mask` types are
`AbstractContainer`'s (Ex: Atom, Residue, Segment, ...)

    ProtoSyn.Mask{T}() where {T <: AbstractContainer}

Creates and empty `Mask` of length 0.

    ProtoSyn.Mask{T}(n::Int) where {T <: AbstractContainer}

Creates a `Mask` of length `n` with all bits set to false.  

Masks can be combined using unary operators, such as `&` and `|` (Only if they
are of the same type. Otherwise, the `promote` function can cast one of the
masks to the other's type, if a corresponding `container` is given).
The `!` unary operator is also available.

*Note:* `ProtoSyn.Mask`s can also de 2D.

# Examples
```jldoctest
julia> ProtoSyn.Mask{Atom}()
ProtoSyn.Mask{Atom}(Bool[])

julia> !ProtoSyn.Mask{Atom}(5)
ProtoSyn.Mask{Atom}(Bool[1, 1, 1, 1, 1])
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
    print(io, "ProtoSyn.Mask{$T}($(size(m.content)))")
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

function tile!(m::Mask{T}) where {T <: AbstractContainer}
    m.content = hcat(collect(Iterators.repeated(m.content, size(m.content)[1]))...)
end

function tile(m::Mask{T}) where {T <: AbstractContainer}
    return Mask{T}(hcat(collect(Iterators.repeated(m.content, size(m.content)[1]))...))
end

Base.adjoint(m::Mask{T}) where {T <: AbstractContainer} = begin
    # ! Note: This is a copy of the original mask.
    # ? Should be view?
    return Mask{T}(m.content')
end

function cross2d(m1::Mask{T}, m2::Mask{T}) where {T <: AbstractContainer}
    _m1 = tile(m1)
    _m2 = tile(m2)
    return _m1 & _m2
end

cross2d(m1::Mask{T}) where {T <: AbstractContainer} = cross2d(m1, m1')