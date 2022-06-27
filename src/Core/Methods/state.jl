# This file should contain functions that work on the system State, such as 
# functions that deal with internal/cartesian coordinate syncs, among others.
using LinearAlgebra

"""
    request_c2i!(state::State; [all::Bool = false])

Sets `state.c2i` to `true`. If `all` is set to `true` (`false`, by default),
update all [`AtomState`](@ref) instances in the given [State](@ref state-types) `state` to
have `:changed` field set to `true`. Return the altered [State](@ref state-types) `state`.

# See also
[`request_i2c!`](@ref) [`c2i!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.request_c2i!(pose.state)
State{Float64}:
 Size: 343
 i2c: false | c2i: true
 Energy: Dict(:Total => Inf)
```
"""
function request_c2i!(state::State; all::Bool = false)
    @info "Requesting cartesian to internal conversion."

    state.c2i = true
    if all
        for atomstate in state.items[4:end]
            atomstate.changed = true
        end
    end
    return state
end


"""
    request_i2c!(state::State; [all::Bool = false])

Sets `state.i2c` to `true`. If `all` is set to `true` (`false`, by default),
update the first [`AtomState`](@ref) instance in the given [State](@ref state-types) `state`
(in the Root) to have `:changed` field set to `true`. Return the altered
[State](@ref state-types) `state`.

# See also
[`request_c2i!`](@ref) [`i2c!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.request_i2c!(pose.state)
State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
```
"""
function request_i2c!(s::State; all::Bool=false)
    @info "Requesting internal to cartesian conversion."
    s[0].changed = all
    s.i2c = true
    s
end


export sync!
"""
    sync!(state::State, topology::Topology)
    
Check whether the given [`State`](@ref) `state` instance has either `i2c` or
`c2i` flag set to `true` and if so update the cartesian/internal coordinates
accordingly. Return the altered [`State`](@ref) instance. 

    sync!(pose::Pose)

Check whether the given [Pose](@ref pose-types) instance has either `i2c` or
`c2i` flag set to `true` in its `pose.state` field and if so update the
cartesian/internal coordinates accordingly. Return the altered [Pose](@ref pose-types) instance. 

!!! ukw "Note:"
    Requesting both `i2c` and `c2i` conversions simultaneously is not possible
    and will result in an error. Consider calling [`i2c!`](@ref) or
    [`c2i!`](@ref) to choose one of the coordinate systems to be synched. 

# See also
[`i2c!`](@ref) [`c2i!`](@ref)

# Examples
```jldoctest
julia> sync!(pose.state, pose.graph)
State{Float64}:
 Size: 343
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
```
"""
function sync!(state::State, topology::Topology)::State

    if state.c2i && state.i2c
        error("unable to request simultaneous i->c and c->i coordinate conversion")
    elseif state.c2i
        @info "Applying pending cartesian to internal changes ..." 
        c2i!(state, topology)
    elseif state.i2c
        @info "Applying pending internal to cartesian changes ..." 
        i2c!(state, topology)
    end

    return state
end

function sync!(pose::Pose)::Pose
    sync!(pose.state, pose.graph)
    return pose
end


"""
    c2i!(state::State{T}, top::Topology)
    
Update the internal coordinates to match the measured cartesian coordinates,
in the given [`State`](@ref) `state`. Note that only the [`AtomState`](@ref)
instances with `:changed` field set to `true` will be updated, and the flag is
therefore changed to `false`. Return the altered [`State`](@ref) `state`
instance. If `state.c2i` is not set to `true`, return the original
[`State`](@ref) `state` instance, without changes.

# See also
[`i2c!`](@ref) [`request_c2i!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.c2i!(pose.state, pose.graph)
State{Float64}:
 Size: 343
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
```
"""
function c2i!(state::State{T}, top::Topology) where T

    state.c2i == false && return state
    for atom in eachatom(top)
        istate = state[atom.id]
        # ? Contrarily to i2c! function, c2i! does not automatically set the
        # ? ascedents :changed field to true. Is this by design or should it be
        # ? altered in future versions?
        !(istate.changed) && continue
        (i, j, k, l) = atom.ascendents
        jstate = state[j]
        kstate = state[k]
        
        # Bond
        istate.b       = ProtoSyn.distance(jstate, istate)

        # Angle
        istate.θ       = ProtoSyn.angle(kstate, jstate, istate)

        # Dihedral
        istate.ϕ       = ProtoSyn.dihedral(state[l], kstate, jstate, istate)
        istate.Δϕ      = 0
        istate.changed = false # * Reset the :changed flag
    end

    state.c2i = false
    return state
end


"""
    i2c!(state::State{T}, top::Topology)
    
Update the cartesian coordinates to match the measured internal coordinates,
in the given [`State`](@ref) `state`. Note that only the [`AtomState`](@ref)
instances with `:changed` field set to `true` will be updated, and the flag is
therefore changed to `false`. Return the altered [`State`](@ref) `state`
instance. If `state.c2i` is not set to `true`, return the original
[`State`](@ref) `state` instance, without changes.

!!! ukw "Note:"
    Any [`AtomState`](@ref) that requires an update (has `:changed` flag set to
    `true`) will cause all downstream residues to be updated as well (in the
    same [Graph](@ref state-types)).

# See also
[`c2i!`](@ref) [`request_i2c!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.i2c!(pose.state, pose.graph)
State{Float64}:
 Size: 343
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
```
"""
function i2c!(state::State{T}, top::Topology) where T
    
    state.i2c == false && return state

    vjk = MVector{3, T}(0, 0, 0)
    vji = MVector{3, T}(0, 0, 0)
    n   = MVector{3, T}(0, 0, 0)
    xi  = MVector{3, T}(0.0, 0.0, 0.0)
    
    queue = Atom[]

    root = ProtoSyn.root(top)
    root_changed = state[root].changed

    for child in root.children
        state[child].changed |= root_changed # Updates state[child].changed to "true" only if 'root_changed' is true.
        push!(queue, child)
    end
    
    while !isempty(queue)
        atom = popfirst!(queue)

        (i, j, k) = atom.ascendents
        istate = state[i]
        
        for child in sort_children(atom)
            # Updates state[child].changed to "true" only if 'istate.changed' is
            # true. (which is, if root_changed is true)
            state[child].changed |= istate.changed
            push!(queue, child)
        end
        !(istate.changed) && continue
        
        jstate = state[j]        
        kstate = state[k]        
        Ri = istate.r
        Rj = jstate.r
        
        # local coord system
        b = istate.b # distance
        sθ, cθ = sincos(istate.θ)  # angle
        sϕ, cϕ = sincos(istate.ϕ + jstate.Δϕ)  # dihedral
        
        x_1 = -b*cθ
        x_2 =  b*cϕ*sθ
        x_3 =  b*sϕ*sθ
        
        # rotate to parent coord system
        @nexprs 3 u -> vji[u] = Rj[u, 1]*x_1 + Rj[u, 2]*x_2 + Rj[u, 3]*x_3
        
        # UPDATE ROTATION MATRIX
        # @nexprs 3 u -> vjk_u = kstate.t[u] - jstate.t[u]
        @. vjk = kstate.t - jstate.t
        
        # column 1 (x)
        @nexprs 3 u -> Ri[u, 1] = vji[u]/b
        
        # column 3 (z)
        @cross u n[u] vji[u] vjk[u]
        dn = sqrt(dot(n,n))
        @nexprs 3 u -> Ri[u, 3] = n[u]/dn
        
        # column 2 (y)
        @cross u Ri[u, 2] Ri[u, 3] Ri[u, 1]
        
        # move to new position
        @. xi          = vji + jstate.t
        istate.t       = xi
        istate.changed = false
    end

    state.i2c = false
    return state
end


export setoffset!
"""
    setoffset!(state::State{T}, at::Atom, default::Number) where {T <: AbstractFloat}

Rotate all sibling dihedrals of [`Atom`](@ref) `atom` in the given
[`State`](@ref) `state` so that the dihedral angle identified by `atom` is equal
to `default`. Set the `i2c` flag to `true` and return the altered
[`State`](@ref) `state`.

# See also
[`LGrammar`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.setoffset!(pose.state, pose.graph[1][2]["CA"], 3.14)
State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
```
"""
function setoffset!(state::State{T}, atom::Atom, default::Number) where {T <: AbstractFloat}

    if hasparent(atom)
        ϕ = state[atom].ϕ - T(default)
        for child in atom.parent.children
            state[child].ϕ -= ϕ
        end
    end

    state.i2c = true
    state
end

export setdihedral!, getdihedral

"""
    setdihedral!(state::State, atom::Atom, value::T) where {T <: AbstractFloat}

Set the dihedral in [`Atom`](@ref) `atom` of [`State`](@ref) `state` to be
exactly `value` (in radians). Automatically requests internal to cartesian
coordinate conversion (by setting `state.i2c` as `true`). Return the altered
[`State`](@ref) `state`.

    setdihedral!(pose::Pose, sele::AbstractSelection, value::T) where {T <: AbstractFloat}

Alternativelly, set the dihedral in the (first) selected [`Atom`](@ref) instance
given by the `AbstractSelection` `sele` in the given [`Pose`](@ref) `pose` to
the provided `value`. Return the altered [`Pose`](@ref) `pose`.

# See also
[`ascendents`](@ref) [`request_i2c!`](@ref) [`getdihedral`](@ref)
[`rotate_dihedral!`](@ref)

# Examples
```jldoctest
julia> setdihedral!(pose.state, pose.graph[1][1][end], Float64(π))
State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
```
"""
@inline setdihedral!(state::State, atom::Atom, value::T) where {T <: AbstractFloat} = begin
    atom2 = atom.ascendents[2]
    state[atom2].Δϕ += value - getdihedral(state, atom)
    ProtoSyn.request_i2c!(state, all = true)
    return state
end

@inline setdihedral!(pose::Pose, sele::AbstractSelection, value::T) where {T <: AbstractFloat} = begin
    atom  = ProtoSyn.promote(sele, Atom)(pose, gather = true)[1]
    atom2 = atom.ascendents[2]
    pose.state[atom2].Δϕ += value - getdihedral(pose.state, atom)
    ProtoSyn.request_i2c!(pose.state, all = true)
    return pose
end


"""
    getdihedral(state::State, atom::Atom)

Get the current dihedral value for [`Atom`](@ref) `atom` of [`State`](@ref)
`state` (in radians, based on the internal coordinates). This value is the sum
of the intrisic dihedral angle `ϕ` and the second ascendent `Δϕ`.

# See also
[`ascendents`](@ref) [`setdihedral!`](@ref) [`dihedral`](@ref)

# Examples
```jldoctest
julia> getdihedral(pose.state, pose.graph[1][2]["N"])
3.141592653589793
```
"""
@inline getdihedral(state::State, atom::Atom) = begin
    atom2 = atom.ascendents[2]
    return state[atom].ϕ + state[atom2].Δϕ
end


export rotate_dihedral!
"""
    rotate_dihedral!(state::State, atom::Atom, value::T) where {T <: AbstractFloat}

Rotate the dihedral in [`Atom`](@ref) `atom` of [`State`](@ref) `state` by
`value` (in radians, __adds__ to the current dihedral angle). Automatically
requests internal to cartesian coordinate conversion (by setting `state.i2c` as
`true`). Return the altered [`State`](@ref) `state`.

# See also
[`ascendents`](@ref) [`request_i2c!`](@ref) [`getdihedral`](@ref)
[`setdihedral!`](@ref)

# Examples
```jldoctest
julia> rotate_dihedral!(pose.state, pose.graph[1][1][end], Float64(π))
State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
```
"""
@inline rotate_dihedral!(state::State, atom::Atom, value::T) where {T <: AbstractFloat} = begin
    atom2 = atom.ascendents[2]
    state[atom2].Δϕ += value
    ProtoSyn.request_i2c!(state, all = true) # ! Unknown bug when all = false
    return state
end


"""
    reindex(state::State)

Re-indexes the whole [`State`](@ref) `state` (excluding the Root), setting
the `:index` field of [`AtomState`](@ref) instances. Return the altered
[`State`](@ref) `state`.

!!! ukw "Note:"
    Since we are altering a field of [`AtomState`](@ref) structs, the `:changed`
    field will automatically be set to `true` and therefore be updated in a
    future [`sync!`](@ref) call (if either `:i2c` or `:c2i` flag in the
    corresponding `state` is set to `true`).

# See also
[`reindex(::Topology, ::Bool)`](@ref)

# Examples
```jldoctest
julia> reindex(pose.state)
State{Float64}:
 Size: 343
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
```
"""
function reindex(state::State)
    for (index, atomstate) in enumerate(state.items[4:end])
        atomstate.index = index
    end
    return state
end