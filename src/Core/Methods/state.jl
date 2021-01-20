# This file should contain functions that work on the system State, such as 
# functions that deal with internal/cartesian coordinate syncs, among others.
using LinearAlgebra

export sync!

"""
    sync!(state::State, top::Topology)
    
Check whether the given `State` instance has either i2c or c2i flag set to true
and update the cartesian/internal coordinates accordingly. Return the altered
`State` instance.


# Examples
```jldoctest
julia> sync(pose.state, pose.graph)
```
"""
function sync!(state::State, top::Topology)::State

    if state.c2i && state.i2c
        error("unable to request simultaneous i->c and c->i coordinate conversion")
    elseif state.c2i
        c2i!(state, top)
    elseif state.i2c
        i2c!(state, top)
    end

    state
end


"""
    c2i!(state::State{T}, top::Topology)
    
Update the internal coordinates to match the measured cartesian coordinates,
in the given `State`. Return the aletered `State` instance.


# Examples
```jldoctest
julia> c2i!(pose.state, pose.graph)
```
"""
function c2i!(state::State{T}, top::Topology) where T

    for atom in eachatom(top)
        (i, j, k, l) = atom.ascendents
        istate = state[i]
        jstate = state[j]
        kstate = state[k]
        
        # Bond
        istate.b = ProtoSyn.distance(jstate, istate)

        # Angle
        istate.θ = ProtoSyn.angle(kstate, jstate, istate)

        # Dihedral
        istate.ϕ = ProtoSyn.dihedral(state[l], kstate, jstate, istate)
        istate.Δϕ = 0
    end

    state.c2i = false
    state
end


"""
    i2c!(state::State{T}, top::Topology)
    
Update the cartesian coordinates to match the current internal coordinates,
in the given `State`. Return the aletered `State` instance.


# Examples
```jldoctest
julia> i2c!(pose.state, pose.graph)
```
"""
function i2c!(state::State{T}, top::Topology) where T
    
    vjk = MVector{3, T}(0, 0, 0)
    vji = MVector{3, T}(0, 0, 0)
    n   = MVector{3, T}(0, 0, 0)
    xi  = MVector{3, T}(0.0, 0.0, 0.0)
    
    queue = Atom[]

    root = origin(top)
    root_changed = state[root].changed

    for child in root.children
        state[child].changed |= root_changed # Updates state[child].changed to "true" only if 'root_changed' is true.
        push!(queue, child)
    end
    
    while !isempty(queue)
        atom = popfirst!(queue)
        (i, j, k) = atom.ascendents
        
        istate = state[i]

        for child in atom.children
            # Updates state[child].changed to "true" only if 'istate.changed' is
            # true. (which is, if root_changed is true)
            state[child].changed |= istate.changed
            push!(queue, child)
        end
        !(istate.changed) && continue
        istate.changed = false
        
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
        @. xi = vji + jstate.t
        # println(istate)
        istate.t = xi
        # println(istate)

        # @. state.x[1:3, i] = xi
    end

    state.i2c = false
    state.c2i = false
    state
end



export setoffset!

"""
    setoffset!(state::State{T}, at::Atom, default::Number) where T

Rotate all sibling dihedrals of `at` so that the dihedral angle identified by
`at` is equal to `default`.
"""
setoffset!(state::State{T}, at::Atom, default::Number) where T = begin

    if hasparent(at)
        ϕ = state[at].ϕ - T(default)
        for child in at.parent.children
            state[child].ϕ -= ϕ
        end
    end
    state.i2c = true
    state
end

export setdihedral!, getdihedral

"""
    setdihedral!(s::State, at::Atom, val::T) where {T <: AbstractFloat}

Set the dihedral in `Atom` `at` of `State` `s` to be `val` (in radians).

# Examples
```jldoctest
julia> setdihedral!(pose.state, pose.graph[1][1][end], π)
```
"""
@inline setdihedral!(s::State, at::Atom, val::T) where {T <: AbstractFloat} = begin
    at2 = at.ascendents[2]
    s[at2].Δϕ += val - getdihedral(s, at)
    ProtoSyn.request_i2c(s)
    s
end


"""
    getdihedral(s::State, at::Atom)

Get the dihedral value for `Atom` `at` of `State` `s` (in radians).
This value is the sum of the intrisic dihedral angle and the second ascendent
Δϕ.

# Examples
```jldoctest
julia> getdihedral(pose.state, pose.graph[1][1][end])
3.1415926535897
```
"""
@inline getdihedral(s::State, at::Atom) = begin
    at2 = at.ascendents[2]
    return s[at].ϕ + s[at2].Δϕ
end

export get_cartesian_matrix
function get_cartesian_matrix(::Type{T}, pose::Pose) where {T <: AbstractFloat}

    matrix = Vector{Vector{T}}()
    for atom_state in pose.state
        push!(matrix, atom_state.t)
    end

    return matrix
end

get_cartesian_matrix(pose::Pose) = get_cartesian_matrix(Units.defaultFloat, pose)


export rotate_dihedral!

"""
    rotate_dihedral!(s::State, at::Atom, val::T) where {T <: AbstractFloat}

Rotate the dihedral in `Atom` `at` of `State` `s` by `val` (in radians).

# Examples
```jldoctest
julia> rotate_dihedral!(pose.state, pose.graph[1][1][end], π)
```
"""
@inline rotate_dihedral!(s::State, at::Atom, val::T) where {T <: AbstractFloat} = begin
    at2 = at.ascendents[2]
    s[at2].Δϕ += val
    ProtoSyn.request_i2c(s)
    s
end

function reindex(s::State)
    for (index, atomstate) in enumerate(s.items[4:end])
        atomstate.index = index
    end
end
