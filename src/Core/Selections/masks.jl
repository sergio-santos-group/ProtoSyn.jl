"""
    ProtoSyn.Mask{T}(content::BitVector) where {T <: AbstractContainer}

A `Mask` is a typed wrapped on a `BitVector`. `Mask` types are `AbstractContainer`'s
(Ex: Atom, Residue, Segment, ...)

    ProtoSyn.Mask{T}() where {T <: AbstractContainer}

Creates and empty `Mask` of length 0.

    ProtoSyn.Mask{T}(n::Int) where {T <: AbstractContainer}

Creates a `Mask` of length `n` with all bits set to false.  

Masks can be combined using unary operators, such as `&` and `|` (Only if they
are of the same type. Otherwise, the `promote` function can cast one of the
masks to the other's type, if a corresponding `container` is given).
The `!` unary operator is also available.

# Examples
```jldoctest
julia> ProtoSyn.Mask{Atom}()
ProtoSyn.Mask{Atom}(Bool[])

julia> !ProtoSyn.Mask{Atom}(5)
ProtoSyn.Mask{Atom}(Bool[1, 1, 1, 1, 1])
```
"""
mutable struct Mask{T <: AbstractContainer}
    content::BitVector
end
Mask{T}() where {T <: AbstractContainer} = Mask{T}(BitVector())
Mask{T}(n::Int) where {T <: AbstractContainer} = Mask{T}(falses(n))

selection_type(::Mask{T}) where {T} = T


# --- Utilities ----------------------------------------------------------------
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