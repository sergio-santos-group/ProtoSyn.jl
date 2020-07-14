mutable struct Mask{T <: AbstractContainer}
    content::BitVector
end
Mask{T}() where {T <: AbstractContainer} = Mask{T}(BitVector())
Mask{T}(n::Int) where {T <: AbstractContainer} = Mask{T}(falses(n))

selection_type(m::Mask{T}) where {T} = T


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