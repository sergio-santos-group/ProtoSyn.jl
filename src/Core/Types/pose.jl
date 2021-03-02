export Pose

"""
    Pose{T <: AbstractContainer}(graph::T, state::State)
    
Return a Pose instance. A Pose is a system structure, having both the
interaction graph and the current state of the system represented.

    Pose(::T, frag::Fragment) where {T <: AbstractContainer}

Return a Pose instance from a fragment, where the State is empty/blank.
"""
mutable struct Pose{T <: AbstractContainer}
    graph::T
    state::State
    Pose(c::T, s::State) where {T <: AbstractContainer}= begin
        c.id != s.id && error("unpairable container (ID: $(c.id)) and state (ID: $(s.id))")
        new{T}(c, s)
    end
end

export Fragment

"""
    Pose{T <: AbstractContainer}(graph::T, state::State)
    
Return a Pose instance. A Pose is a system structure, having both the
interaction graph and the current state of the system represented.

    Pose(::T, frag::Fragment) where {T <: AbstractContainer}

Return a Pose instance from a Fragment, where the State is empty/blank.

!!! note
    A Fragment is a Pose{Segment}.

*See also:*
    
`fragment` `Builder.fragment`

"""
const Fragment = Pose{Segment}

Pose(::Type{T}, frag::Fragment) where {T <: AbstractFloat} = begin
    frag2 = copy(frag)
    top = Topology(frag2.graph.name, 1)
    state = State{T}()
    state.id = top.id
    pose = Pose(top, state)
    Base.append!(pose, frag2)

    ProtoSyn.request_i2c(state; all=true)
    return pose
end

Pose(frag::Fragment) = Pose(Float64, frag)


export fragment

"""
    fragment(pose::Pose{Topology})
    
Return a fragment from a given Pose `pose`. The pose must have a single segment.

# Examples
```jldoctest
julia> frag = fragment(pose)
```
"""
function fragment(pose::Pose{Topology})
    
    length(pose.graph) != 1 && error("only topologies with a single segment can be turned into fragments")
    
    topology = pose.graph
    segment = topology[1]
    state = splice!(pose.state, 1:count_atoms(segment))
    detach(segment)
    segment.id = state.id = genid()
    segment.name = topology.name
    segment.container = nothing

    Pose(segment, state)
end