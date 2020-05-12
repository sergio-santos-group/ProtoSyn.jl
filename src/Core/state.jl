using LinearAlgebra: I
using StaticArrays

#region AtomState --------------------------------------------------------------

export AtomState
mutable struct AtomState{T<:AbstractFloat}
    t::MVector{3,T}         # translation vector
    r::MMatrix{3,3,T,9}     # local to global rotation matrix
    
    # internal coordinates
    b::T                    # self<->parent bond length
    θ::T                    # self<->parent<->grandparent angle
    ϕ::T                    # self<->parent<->grandparent<->grand-grandparent dihedral
    
    Δϕ::T                   # change in dihedral angles (to be applied to children)
    changed::Bool           # flag
end

AtomState{T}() where {T} = begin
    AtomState{T}(
        MVector{3,T}(zeros(T, 3)),
        MMatrix{3,3,T,9}(Matrix{T}(I,3,3)),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        false
    )
end

Base.setproperty!(ns::AtomState{T}, key::Symbol, val) where T = begin
    if key == :Δϕ
        setfield!(ns, :changed, true)
    end
    setfield!(ns, key, T(val))
end


#endregion AtomState


#region State ------------------------------------------------------------------


export State
mutable struct State{T<:AbstractFloat}
    items::Vector{AtomState{T}} 
    size::Int           # number of items
    id::Int             # id to be matched to the corresponding container
    i2c::Bool           # flag to request internal to cartesian conversion
    c2i::Bool           # flag to request cartesian to internal conversion
    index_offset::Int   # ...

    # x::Array{T,2}
    # f::Array{T,2}
    # v::Array{T,2}
end

State{T}(n::Int) where T = begin
    items = Vector{AtomState{T}}(undef, n+3)
    for i=1:n+3
        items[i] = AtomState{T}()
    end
    items[1].t[1] = -1.0
    items[1].t[2] =  1.0
    items[2].t[1] = -1.0
    State{T}(items, n, -1, false, false, 3)
end

State(n::Int) = State{Float64}(n)

State(::Type{T}, n::Int) where T = State{T}(n)

State(items::Vector{AtomState{T}}) where T = begin
    s = State{T}(items, length(items), -1, false, false, 0)
end
State{T}() where T = State{T}(0)

Base.getindex(s::State, i::Int) = begin
    s.items[i+s.index_offset]
end
Base.getindex(s::State, at::Atom) = begin
    s.items[at.index+s.index_offset]
end

Base.firstindex(s::State) = 1-s.index_offset
Base.lastindex(s::State) = s.size

Base.eltype(::Type{State{T}}) where T = T

Base.splice!(s::State, range::UnitRange{Int}) = begin
    i = s.index_offset + range.start
    j = s.index_offset + range.stop
    s2 = State(splice!(s.items, i:j))
    s.size -= s2.size
    s2
end

Base.append!(s1::State{T}, s2::State{T}) where T = begin
    for i = 1:s2.size
        push!(s1.items, s2[i])
    end
    s1.size += s2.size
    s1
end

Base.insert!(s1::State{T}, index::Integer, s2::State{T}) where T = begin
    index += s1.index_offset
    for i = 0:s2.size-1
        insert!(s1.items, index+i, s2[i+1])
    end
    s1.size += s2.size
    s1
end


Base.copy(s::State{T}) where T = deepcopy(s)


request_c2i(s::State; all=false) = (s.c2i = true)

request_i2c(s::State; all=false) = begin
    s[0].changed = all
    s.i2c = true
    s
end


#endregion State