using LinearAlgebra: I
using StaticArrays
using Printf

#region AtomState --------------------------------------------------------------

export AtomState, StateMatrix

"""
    AtomState{T}(parent::Any, index::Int, t::MVector{3, T}, r::MMatrix{3, 3, T, 9}, b::T, θ::T, ϕ::T, δ::T, Δϕ::T, changed::Bool) where {T <: AbstractFloat}
    
An [`AtomState`](@ref) instance. Holds information regarding the state of the atom,
including the cartesian and internal coordinates.

    AtomState([::T]) where {T <: AbstractFloat}

Return an empty [`AtomState`](@ref) instance, with all default values.

!!! ukw "Note:"
    If no type `T <: AbstractFloat` is provided, `Units.defaultFloat` will be used.

!!! ukw "Note:"
    Changing this [`AtomState`](@ref) cartesian coordinates will update and sync
    changes with a sibling [`StateMatrix`](@ref), as long as this structs `parent`
    is correctly set.

# Fields
* `parent::Any` - If set, should point to the parent [`State`](@ref) containing this [`AtomState`](@ref) (default: nothing)
* `index::Int` - the index of this [`AtomState`](@ref)
* `t::MVector{3, T}` - The translation vector from origin (cartesian coordinates) (default: empty)
* `r::MMatrix{3, 3, T, 9}` - The rotation matrix for cartesian to internal coordinate sync (default: empty)
* `b::T` - Distance (in Angstrom Å) to parent atom (default: 0)
* `θ::T` - Angle (in radians) to ascendent atoms (default: 0)
* `ϕ::T` - Dihedral angle (in radians) to ascendent atoms (default: 0)
* `δ::T` - Atomic partial charge (default: 0.0)
* `Δϕ::T` - Dihedral angle change (in radians) to be applied to children atoms (default: 0)
* `changed::Bool` - Flag indicating whether this [`AtomState`](@ref) has been modified (useful in some functions such as [`i2c!`](@ref ProtoSyn.i2c!) and [`c2i!`](@ref ProtoSyn.c2i!)) (default: false)

!!! ukw "Note:"
    The `Δϕ` field in [`AtomState`](@ref) allows for easy set-up of dihedral angles in molecular structures.
    By modifying `Δϕ` in an atom, all child atoms will be rotated by the same amount, even in being branched structures.

# See also
[`State`](@ref) [`StateMatrix`](@ref)

# Examples
```jldoctest
julia> AtomState()
AtomState{Float64}:
 Index: -1
 T: [0.000, 0.000, 0.000]
 b: 0.000 Å | θ:  0.000 rad (   0.00°) | ϕ:  0.000 rad (   0.00°) | Δϕ:  0.000 rad (   0.00°)
 δ: 0.000
 Changed: false
```
"""
mutable struct AtomState{T <: AbstractFloat}
    parent
    index::Int
    t::MVector{3, T}        # translation vector
    r::MMatrix{3, 3, T, 9}  # local to global rotation matrix
    
    # internal coordinates
    b::T                    # self<->parent bond length
    θ::T                    # self<->parent<->grandparent angle
    ϕ::T                    # self<->parent<->grandparent<->grand-grandparent dihedral

    # charge
    δ::T                    # partial charge
    
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
        if ns.parent !== nothing
            update_state_matrix(ns.parent.x, val, :, ns.index, update_items = false)
        end
    else
        setfield!(ns, :changed, true)
        setfield!(ns, key, val)
    end
    ns
end

Base.findfirst(x::AtomState{T}, v::Vector{AtomState{T}}) where {T <: AbstractFloat} = begin
    return findfirst(i -> i === x, v)
end

function Base.show(io::IO, as::AtomState{T}) where {T <: AbstractFloat}
    println(io, "AtomState{$T}:")
    println(io, " Index: $(as.index)")
    @printf(io, " T: [%.3f, %.3f, %.3f]\n", as.t[1], as.t[2], as.t[3])
    # @printf(io, " b: $(as.b) | θ: $(as.θ) rad ($(deg2rad(as.θ))°) | ϕ: $(as.ϕ) | Δϕ: $(as.Δϕ)")
    @printf(io, " b: %5.3f Å | θ: %6.3f rad (%7.2f°) | ϕ: %6.3f rad (%7.2f°) | Δϕ: %6.3f rad (%7.2f°)\n", as.b, as.θ, rad2deg(as.θ), as.ϕ, rad2deg(as.ϕ), as.Δϕ, rad2deg(as.Δϕ))
    @printf(io, " δ: %5.3f\n", as.δ)
    println(io, " Changed: $(as.changed)")
end

#endregion AtomState

"""
    StateMatrix{T}(parent::Any, coords::Matrix{T}) where {T <: AbstractFloat}

A [`StateMatrix`](@ref) instance holds the cartesian coordinates of all
[`AtomState`](@ref) instances in a [`State`](@ref). It is an overload of a
`Matrix` struct, and used to efficiently apply substantial changes to a system
or a large number of atoms. 

!!! ukw "Note:"
    If no type `T <: AbstractFloat` is provided, `Units.defaultFloat` will be
    used.

!!! ukw "Note:"
    Changing this [`StateMatrix`](@ref) cartesian coordinates will update and
    sync changes with a sibling [`AtomState`](@ref), as long as this structs
    `parent` is correctly set. In such case, the `parent.c2i` flag is also set
    to `true`.

# Fields
* `parent::Any` - If set, should point to the parent [`State`](@ref) containing this [`StateMatrix`](@ref) (default: nothing)
* `coords::Matrix` - The cartesian coordinates matrix

# See also
[`State`](@ref) [`AtomState`](@ref)

# Examples
```jldoctest
julia> StateMatrix(zeros(3, 3))
StateMatrix{Float64}:
 Parent set: false
3×3 Matrix{Float64}:
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0
```
"""
mutable struct StateMatrix{T <: AbstractFloat}
    parent
    coords::Matrix{T}
end

StateMatrix(x::Matrix{T}) where {T <: AbstractFloat} = StateMatrix{T}(nothing, x)

Base.length(x::StateMatrix{T}) where {T <: AbstractFloat} = length(x.coords)
Base.size(x::StateMatrix{T}) where {T <: AbstractFloat} = size(x.coords)

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
    Base.show(io, MIME"text/plain"(), sm.coords)
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

Base.lastindex(x::StateMatrix{T}, i::Int64) where {T <: AbstractFloat} = begin
    return size(x.coords)[2]
end

Base.iterate(x::StateMatrix{T}, (a,)=(1,)) where {T <: AbstractFloat} = begin
    if a > size(x)[2]
        return nothing
    end
    (x.coords[:, a], (a+1,))
end

#region State ------------------------------------------------------------------

export State

"""
    State([::T]) where {T <: AbstractFloat}

A `State` is a complete representation of a given molecular system regarding the
position, the forces felt and the resulting energetic contribution by each atom.
Each [`AtomState`](@ref) in a system.

    State([::T], n::Int) where {T <: AbstractFloat}

Return a [`State`](@ref) with size `n`, with all `items` set to be an empty
[`AtomState`](@ref).

    State([::T], items::Vector{AtomState{T}}) where {T <: AbstractFloat}

Return a [`State`](@ref) with size `length(items)`, with all the given `items`.

!!! ukw "Note:"
    If no type `T <: AbstractFloat` is provided, `Units.defaultFloat` will be used.

# Fields
* `items::Vector{AtomState{T}}` - The list of [`AtomState`](@ref) instances in this [`State`](@ref) (default: empty)
* `size::Int` - The number of [`AtomState`](@ref) instances in this [`State`](@ref) (default: 0)
* `id::Int` - The Id (to be matched to the corresponding [`Topology`](@ref) in a [`Pose`](@ref)) (default: -1)
* `i2c::Bool` - Flag indicating this [`State`](@ref) needs to be synched to cartesian coordinates (default: false)
* `c2i::Bool` - Flag indicating this [`State`](@ref) needs to be synched to internal coordinates (default: false)
* `index_offset::Int` - Ignore the first *N* [`AtomState`](@ref) instances (default: 3)
* `x::StateMatrix{T}` - The cartesian coordinate matrix of all [`AtomState`](@ref) instances in this [`State`](@ref) (default: empty)
* `f::Matrix{T}` - The force matrix of all [`AtomState`](@ref) instances in this [`State`](@ref) (default: empty)
* `e::Dict{Symbol, T}` - The list of all energetic components evaluated for this [`State`](@ref) (default: empty)

# See also
[`Topology`](@ref) [`AtomState`](@ref) [`StateMatrix`](@ref)

# Examples
```jldoctest
julia> State(4)
State{Float64}:
 Size: 4
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
```
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

    state = State{T}(items, n, -1, false, false, 3, StateMatrix(zeros(T, 3, n)), zeros(T, 3, n), Dict(:Total => Inf))
    for (index, item) in enumerate(state)
        item.parent = state
        item.index  = index
    end
    state.x.parent = state

    return state
end

State(n::Int) = State{Units.defaultFloat}(n)
State(::Type{T}, n::Int) where T = State{T}(n)
State(::Type{T}, items::Vector{AtomState{T}}) where {T <: AbstractFloat} = begin
    state = State{T}(items, length(items), -1, false, false, 0, StateMatrix(zeros(T, 3, length(items))), zeros(T, 3, length(items)), Dict(:Total => Inf))
    for (index, item) in enumerate(state.items)
        item.parent = state
        item.index  = index
        item.t      = item.t # Force update of state.x
    end
    state.x.parent = state

    return state
end

State(items::Vector{AtomState{T}}) where T = State(T, items)

State(::Type{T}) where {T <: AbstractFloat} = State{T}(0)
State() = State{Units.defaultFloat}(0)
State{T}() where {T <: AbstractFloat} = State{T}(0)

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

Base.getindex(s::State, ats::Vector{Atom}) = begin
    tr = Vector{AtomState}()
    for at in ats
        push!(tr, s[at])
    end
    
    return tr
end

Base.getindex(s::State, seg::Segment) = begin
    s.items[(seg[1][1].index + s.index_offset):(seg[end][end].index + s.index_offset)]
end

function Base.setindex!(s::State, ns::AtomState{T}, at::Atom) where {T <: AbstractFloat}
    os         = s[at]
    os.parent  = s
    os.b       = ns.b
    os.θ       = ns.θ
    os.ϕ       = ns.ϕ
    os.Δϕ      = ns.Δϕ
    os.δ       = ns.δ
    os.changed = ns.changed
    os.index   = ns.index
    os.t       = ns.t
    os.r       = ns.r

    return os
end

Base.firstindex(s::State) = 1-s.index_offset
Base.lastindex(s::State) = s.size

Base.eltype(::Type{State{T}}) where T = T
Base.eltype(::State{T}) where T = T

Base.splice!(s::State{T}, range::UnitRange{Int}) where {T <: AbstractFloat}= begin
    # Both range ends are included in the spliced region
    i = s.index_offset + range.start
    j = s.index_offset + range.stop
    s2 = State(T, splice!(s.items, i:j))
    downstream = s.x.coords[:, 1:(range.start-1)]
    upstream   = s.x.coords[:, (range.stop+1):end]
    s.x.coords = hcat(downstream, upstream)
    s.size -= s2.size
    s2
end

Base.splice!(s::State, index::Int) = begin
    _index     = s.index_offset + index
    s2         = State([splice!(s.items, _index)]) # * Already sets state.x
    s2.f      .= s.f[:, index]
    s.x.coords = hcat(s.x.coords[:, 1:(index-1)], s.x.coords[:, (index+1):end])
    s.f        = hcat(s.f[:, 1:(index-1)], s.f[:, (index+1):end])
    s.size    -= s2.size
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
    s1.x.coords = hcat(s1.x.coords, s2.x.coords)
    s1
end

Base.insert!(s1::State{T}, index::Integer, s2::State{T}) where T = begin
    # * Note: state.x.coords doesn't include the 3 origin atoms.
    s1.x.coords = hcat(s1.x.coords[:, 1:(index-1)], s2.x.coords, s1.x.coords[:, index:end])
    # * Note: state.f doesn't include the 3 origin atoms.
    s1.f = hcat(s1.f[:, 1:(index-1)], s2.f, s1.f[:, index:end])
    # * Now we can add the index_offset
    index += s1.index_offset
    for i = 0:s2.size-1
        insert!(s1.items, index+i, s2[i+1])
        # * The items need to have the parent state re-assigned
        s1.items[index+i].parent = s1
    end
    s1.size += s2.size
    s1
end

Base.copy(as::AtomState{T}) where {T <: AbstractFloat} = begin
    return ProtoSyn.AtomState(
        nothing, as.index, copy(as.t), copy(as.r),
        as.b, as.θ, as.ϕ, as.δ, as.Δϕ, as.changed)
end

Base.push!(s::State{T}, as::AtomState{T}) where T = begin
    push!(s.items, as)
    s.x.coords = hcat(s.x.coords, as.t)
    s.f = hcat(s.f, [T(0.0), T(0.0), T(0.0)])
    s.size += 1
    as.parent = s
    return s
end

#endregion State


function atmax(state::State{T}, comp::Symbol) where T
    m   = T(0)
    idx = 0
    x   = getproperty(state, comp)
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