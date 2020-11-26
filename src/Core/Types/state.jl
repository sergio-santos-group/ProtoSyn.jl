using LinearAlgebra: I
using StaticArrays

#region AtomState --------------------------------------------------------------

export AtomState

"""
    AtomState{T <: AbstractFloat}()

An `AtomState` instance. Holds information regarding the state of the atom,
including the cartesian and internal coordinates.
"""
mutable struct AtomState{T <: AbstractFloat}
    parent
    index::Int
    t::MVector{3,T}         # translation vector
    r::MMatrix{3,3,T,9}     # local to global rotation matrix
    
    # internal coordinates
    b::T                    # self<->parent bond length
    θ::T                    # self<->parent<->grandparent angle
    ϕ::T                    # self<->parent<->grandparent<->grand-grandparent dihedral <- GRANDCHILDREN ?!
    
    Δϕ::T                   # change in dihedral angles (to be applied to children)
    changed::Bool           # flag
end

AtomState{T}() where {T} = begin
    AtomState{T}(
        nothing,
        -1,
        MVector{3,T}(zeros(T, 3)),
        MMatrix{3,3,T,9}(Matrix{T}(I,3,3)),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        false)
end

AtomState() = AtomState{Units.defaultFloat}()

Base.setproperty!(ns::AtomState{T}, key::Symbol, val) where T = begin
    # Intercepts `setproperty!` to change the :changed field to true.

    if key == :Δϕ
        setfield!(ns, :changed, true)
        setfield!(ns, key, T(val))
    elseif key == :changed
        setfield!(ns, :changed, val)
    elseif key == :t
        setfield!(ns, :changed, true)
        setfield!(ns, key, MVector{3, T}(val))
        update_state_matrix(ns.parent.x, val, :, ns.index, update_items = false)
        ns.parent.c2i = true
    else
        setfield!(ns, :changed, true)
        setfield!(ns, key, val)
    end
    ns
end

function Base.show(io::IO, as::AtomState{T}) where {T <: AbstractFloat}
    println(io, "AtomState{$T}:")
    println(io, " Index: $(as.index)")
    println(io, " T: $(as.t)")
    println(io, " b: $(as.b) | θ: $(as.θ) | ϕ: $(as.ϕ) | Δϕ: $(as.Δϕ)")
    println(io, " Changed: $(as.changed)")
end

#endregion AtomState


mutable struct StateMatrix{T <: AbstractFloat}
    parent
    coords::Matrix{T}
end

StateMatrix(x::Matrix{T}) where {T <: AbstractFloat} = StateMatrix{T}(nothing, x)

Base.getindex(x::StateMatrix{T}, i::Int) where {T <: AbstractFloat} = x.coords[i]
Base.getindex(x::StateMatrix{T}, i::UnitRange{<: Real}) where {T <: AbstractFloat} = x.coords[i]
Base.getindex(x::StateMatrix{T}, i::UnitRange{<: Real}, j::Int) where {T <: AbstractFloat} = x.coords[i, j]
Base.getindex(x::StateMatrix{T}, j::Int, i::Colon) where {T <: AbstractFloat} = x.coords[j, i]
Base.getindex(x::StateMatrix{T}, i::Colon, j::Int) where {T <: AbstractFloat} = x.coords[i, j]
Base.getindex(x::StateMatrix{T}, i::Colon, j::UnitRange{<: Real}) where {T <: AbstractFloat} = x.coords[i, j]
Base.getindex(x::ProtoSyn.StateMatrix{T}, i::Colon) where {T <: AbstractFloat} = x.coords[i]
Base.getindex(x::ProtoSyn.StateMatrix{T}, i::Colon, j::Array{Int64,1}) where {T <: AbstractFloat} = x.coords[i, j]

function Base.show(io::IO, sm::StateMatrix{T}) where {T <: AbstractFloat}
    println(io, "StateMatrix{$T}:")
    println(io, " Parent set: $(sm.parent !== nothing)")
    show(io, MIME"text/plain"(), sm.coords)
end

function update_state_matrix(x::StateMatrix{T}, val, i::Colon, j::Int; update_items = true) where {T <: AbstractFloat}
    setindex!(x.coords, val, i, j)
    if update_items && x.parent !== nothing
        x.parent[j].t = MVector{3, T}(val)
    end
end

function update_state_matrix(x::StateMatrix{T}, vals, i::Colon, j::Union{UnitRange{<: Real}, Array{<: Real, 1}}; update_items = true) where {T <: AbstractFloat}
    setindex!(x.coords, vals, i, j)
    if update_items && x.parent !== nothing
        for (value_index, state_index) in enumerate(j)
            x.parent[state_index].t = MVector{3, T}(vals[i, value_index])
        end
    end
end

Base.setindex!(x::StateMatrix{T}, val::Vector{T}, i::Colon, j::Int) where {T <: AbstractFloat} = begin
    update_state_matrix(x, val, i, j)
end

Base.setindex!(x::StateMatrix{T}, vals::Matrix{T}, i::Colon, j::Union{UnitRange{<: Real}, Array{<: Real,1}}) where {T <: AbstractFloat} = begin
    update_state_matrix(x, vals, i, j)
end



#region State ------------------------------------------------------------------

export State

"""
    State{T <: AbstractFloat}()

A `State` includes all information of all `AtomStates` in a system.

    State{T <: AbstractFloat}(n::Int)

Returns a `State` with size `n`, with all `items` set to be an empty
`AtomState`.

State{T <: AbstractFloat}(items::Vector{AtomState{T}})

Returns a `State` with size `length(items)`, with all the given `items`.
"""
mutable struct State{T <: AbstractFloat}
    items::Vector{AtomState{T}} 
    size::Int           # number of items
    id::Int             # id to be matched to the corresponding container
    i2c::Bool           # flag to request internal to cartesian conversion
    c2i::Bool           # flag to request cartesian to internal conversion
    index_offset::Int   # serves to offset the 3 initial root atoms

    x::StateMatrix{T}
    f::Matrix{T}
    e::Dict{Symbol, T}
end

function Base.show(io::IO, state::State{T}) where {T <: AbstractFloat}
    println(io, "State{$T}:")
    println(io, " Size: $(state.size)")
    println(io, " i2c: $(state.i2c) | c2i: $(state.c2i)")
    println(io, " Energy: $(state.e)")
end

State{T}(n::Int) where T = begin

    items = Vector{AtomState{T}}(undef, n+3)

    for i = 1:n+3
        items[i] = AtomState{T}()
    end

    items[1].t[1] = -1.0
    items[1].t[2] =  1.0
    items[2].t[1] = -1.0

    state = State{T}(items, n, -1, false, false, 3, StateMatrix(zeros(T, 3, n)), zeros(T, 3, n), Dict())
    for (index, item) in enumerate(state)
        item.parent = state
        item.index  = index
    end
    state.x.parent = state

    return state
end

State(n::Int) = State{Units.defaultFloat}(n)
State(::Type{T}, n::Int) where T = State{T}(n)
State(::Type{T}, items::Vector{AtomState{T}}) where T = begin
    state = State{T}(items, length(items), -1, false, false, 0, StateMatrix(zeros(T, 3, length(items))), zeros(T, 3, length(items)), Dict())
    for (index, item) in enumerate(state.items)
        item.parent = state
        item.index  = index
        item.t      = item.t # Force update of state.x
    end
    state.x.parent = state

    return state
end

State(items::Vector{AtomState{T}}) where T = State(T, items)

State{T}() where T = State{T}(0)

# ---

Base.length(st::State{T}) where {T <: AbstractFloat} = st.size

Base.iterate(st::State{T}, (s,)=(1,)) where {T <: AbstractFloat} = begin
    if s > length(st)
        return nothing
    end
    (st.items[s + st.index_offset], (s+1,))
end

Base.getindex(s::State, i::Int) = begin
    s.items[i + s.index_offset]
end
Base.getindex(s::State, at::Atom) = begin
    s.items[at.index + s.index_offset]
end
Base.getindex(s::State, seg::Segment) = begin
    s.items[(seg[1][1].index + s.index_offset):(seg[end][end].index + s.index_offset)]
end

Base.firstindex(s::State) = 1-s.index_offset
Base.lastindex(s::State) = s.size

Base.eltype(::Type{State{T}}) where T = T
Base.eltype(::State{T}) where T = T

Base.splice!(s::State{T}, range::UnitRange{Int}) where {T <: AbstractFloat}= begin
    i = s.index_offset + range.start
    j = s.index_offset + range.stop
    s2 = State(T, splice!(s.items, i:j))
    s.size -= s2.size
    s2
end

Base.splice!(s::State, index::Int) = begin
    _index = s.index_offset + index
    s2 = State([splice!(s.items, _index)])
    s.size -= s2.size
    s2
end

Base.append!(s1::State{T}, s2::State{T}) where T = begin
    for i = 1:s2.size
        push!(s1.items, s2[i])
        s1.items[end].parent = s1
        s1.items[end].index += s1.size
    end
    s1.size += s2.size
    s1.f = hcat(s1.f, zeros(T, 3, s2.size))
    s1.x.coords = hcat(s1.x.coords, zeros(T, 3, s2.size))
    s1
end

Base.insert!(s1::State{T}, index::Integer, s2::State{T}) where T = begin
    index += s1.index_offset
    for i = 0:s2.size-1
        insert!(s1.items, index+i, s2[i+1])
    end
    s1.x.coords = hcat(s1.x.coords[:, 1:(index-1)], s2.x.coords, s1.x.coords[:, index:end])
    # TO DO :: FIX STATE.FORCES (INSERT ZEROS)
    s1.size += s2.size
    s1
end

request_c2i(s::State; all=false) = (s.c2i = true; s)

request_i2c(s::State; all::Bool=false) = begin
    s[0].changed = all
    s.i2c = true
    s
end

#endregion State


#--------------------------------- ???
# abstract type AbstractRepresentation end
# struct CartesianRepresentation{T} <: AbstractRepresentation
#     state::State{T}
# end
# struct InternalRepresentation{T} <: AbstractRepresentation
#     state::State{T}
# end

# export save, save!, restore!, AbstractRepresentation, CartesianRepresentation, InternalRepresentation

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


# function save(x::InternalRepresentation{T}) where T
#     m = Matrix{T}(undef, x.state.size, 3)
#     save!(m, x)
# end

# function save!(m::Matrix{T}, x::InternalRepresentation{T}) where T
#     state = x.state
#     for i=1:state.size
#         s = state[i]
#         m[i, 1] = s.b
#         m[i, 2] = s.θ
#         m[i, 3] = s.ϕ
#     end
#     m
# end

# function restore!(x::InternalRepresentation{T}, m::Matrix{T}) where T
#     state = x.state
#     for i=1:state.size
#         s = state[i]
#         s.b = m[i, 1]
#         s.θ = m[i, 2]
#         s.ϕ = m[i, 3]
#     end
#     x
# end

# export prepare
# function prepare(x::AbstractRepresentation)
#     state = x.state
#     # state.x = zeros(T, state.size, 3)
#     state.f = zeros(eltype(state), state.size, 3)
# end

# Base.reset(s::State{T}) where T = begin
#     e = s.e
#     for k in keys(e); e[k] = 0; end
#     fill!(s.f, 0)
#     s
# end

# export energy,setenergy!
# energy(s::State) = sum(values(s.e))
# energy(s::AbstractRepresentation) = energy(s.state)
# setenergy!(s::State, key::Symbol, value) = (s.e[key]=value; s)


function atmax(state::State{T}, comp::Symbol) where T
    m = T(0)
    idx = 0
    x = getproperty(state, comp)
    for i = 1:state.size
        f = x[1, i]^2 + x[2, i]^2 + x[3, i]^2
        if f > m
            m = f
            idx = i
        end
    end
    return sqrt(m), idx
end

function reset_forces!(state::State{T}) where {T <: AbstractFloat}
    for i in 1:state.size
        for j in 1:3
            state.f[j, i] = T(0)
        end
    end
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