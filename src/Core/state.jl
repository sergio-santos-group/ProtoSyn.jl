using LinearAlgebra: I
using StaticArrays

#region AtomState --------------------------------------------------------------

export AtomState

"""
    AtomState{T<:AbstractFloat}


"""
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
        setfield!(ns, key, T(val))
    elseif key == :changed
        setfield!(ns, :changed, val)
    else
        setfield!(ns, key, T(val))
    end
    ns
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

    x::Opt{Matrix{T}}
    f::Opt{Matrix{T}}
    # f::Array{T,2}
    # v::Array{T,2}
    e::Dict{Symbol,T}
end

State{T}(n::Int) where T = begin
    items = Vector{AtomState{T}}(undef, n+3)
    for i=1:n+3
        items[i] = AtomState{T}()
    end
    items[1].t[1] = -1.0
    items[1].t[2] =  1.0
    items[2].t[1] = -1.0
    # State{T}(items, n, -1, false, false, 3)
    State{T}(items, n, -1, false, false, 3, nothing, nothing, Dict())
end

State(n::Int) = State{Float64}(n)

State(::Type{T}, n::Int) where T = State{T}(n)

State(items::Vector{AtomState{T}}) where T = begin
    s = State{T}(items, length(items), -1, false, false, 0, nothing, nothing, Dict())
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
Base.eltype(::State{T}) where T = T

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


request_c2i(s::State; all=false) = (s.c2i = true; s)

# QUESTION
# What is the "all" keyword for? s[0] is the root?
request_i2c(s::State; all::Bool=false) = begin
    s[0].changed = all
    s.i2c = true
    s
end

# function collectx(s::State{T}) where T
#     x = Matrix{T}(undef, s.size, 3)
#     collectx!(x, s)
# end

# function collectx!(x::Matrix{T}, s::State{T}) where T
#     for i=1:s.size
#         @. x[i,:] = s[i].t
#     end
#     x
# end

#endregion State


#---------------------------------
abstract type AbstractRepresentation end
struct CartesianRepresentation{T} <: AbstractRepresentation
    state::State{T}
end
struct InternalRepresentation{T} <: AbstractRepresentation
    state::State{T}
end

export save, save!, restore!, AbstractRepresentation, CartesianRepresentation, InternalRepresentation

# function backup(x::CartesianView{T}) where T
#     state = x.state
#     m = Matrix{T}(undef, state.size, 3)
#     for i=1:state.size
#         t = state[i].t
#         m[i,1] = t[1]
#         m[i,2] = t[2]
#         m[i,3] = t[3]
#     end
#     m
# end

# function backup!(m::Matrix{T}, x::CartesianView{T}) where T
#     state = x.state
#     for i=1:state.size
#         t = state[i].t
#         m[i,1] = t[1]
#         m[i,2] = t[2]
#         m[i,3] = t[3]
#     end
#     m
# end

# function backup!(x::CartesianView{T}, m::Matrix{T}) where T
#     state = x.state
#     for i=1:state.size
#         t = state[i].t
#         t[1] = m[i,1]
#         t[2] = m[i,2]
#         t[3] = m[i,3]
#     end
#     x
# end


function save(x::InternalRepresentation{T}) where T
    m = Matrix{T}(undef, x.state.size, 3)
    save!(m, x)
end

function save!(m::Matrix{T}, x::InternalRepresentation{T}) where T
    state = x.state
    for i=1:state.size
        s = state[i]
        m[i,1] = s.b
        m[i,2] = s.θ
        m[i,3] = s.ϕ
    end
    m
end

function restore!(x::InternalRepresentation{T}, m::Matrix{T}) where T
    state = x.state
    for i=1:state.size
        s = state[i]
        s.b = m[i,1]
        s.θ = m[i,2]
        s.ϕ = m[i,3]
    end
    x
end

export prepare
function prepare(x::AbstractRepresentation)
    state = x.state
    T = eltype(state)
    state.x = zeros(T, state.size, 3)
    state.f = zeros(T, state.size, 3)
end

Base.reset(s::State{T}) where T = begin
    e = s.e
    for k in keys(e); e[k] = 0; end
    fill!(s.f, 0)
    s
end

export energy,setenergy!
energy(s::State) = sum(values(s.e))
energy(s::AbstractRepresentation) = energy(s.state)
setenergy!(s::State, key::Symbol, value) = (s.e[key]=value; s)


function atmax(state::State{T}, comp::Symbol) where T
    m = T(0)
    idx = 0
    x = getproperty(state, comp)
    for i = 1:state.size
        f = x[i,1]^2 + x[i,2]^2 + x[i,3]^2
        if f > m
            m = f
            idx = i
        end
    end
    return sqrt(m),idx
end

struct NamedView{T}
    parent::State{T}
    index::Int
end

Base.getproperty(nv::NamedView, key::Symbol) = begin
    if key === :b
        return nv.parent.i[nv.index,1]
    elseif key === :θ
        return nv.parent.i[nv.index,2]
    elseif key === :ϕ
        return nv.parent.i[nv.index,3]
    elseif key === :x
        return nv.parent.x[nv.index,1]
    elseif key === :y
        return nv.parent.x[nv.index,2]
    elseif key === :z
        return nv.parent.x[nv.index,3]
    elseif key === :changed
        return nv.parent.changed[nv.index]
    end
    getfield(nv, key)
end

# Base.setproperty!(ns::AtomState{T}, key::Symbol, val) where T = begin
#     if key == :Δϕ
#         setfield!(ns, :changed, true)
#         setfield!(ns, key, T(val))
#     elseif key == :changed
#         setfield!(ns, :changed, val)
#     else
#         setfield!(ns, key, T(val))
#     end
#     ns
# end
