export Pose

"""
    Pose{T <: AbstractContainer}(graph::T, state::State)
    
Return a Pose instance. A Pose is a complete description of a molecular system 
at any given point, having both the interaction graph and the current state of
the system coordinates represented. A [`Pose`](@ref) is typed by an
[AbstractContainer], usually a [`Topology`](@ref).

    Pose(::T, frag::Fragment) where {T <: AbstractContainer}

Return a [Pose](@ref pose-types) instance from a [`Fragment`](@ref), where the [State](@ref state-types) is empty/blank. The
graph contents are copied to the new [Pose](@ref pose-types). If no type
`T <: AbstractFloat` is provided, the `Units.defaultFloat` will be used.

# See also
[`Fragment`](@ref)
"""
mutable struct Pose{T <: AbstractContainer}
    graph::T
    state::State
    Pose(c::T, s::State) where {T <: AbstractContainer} = begin
        c.id != s.id && error("unpairable container (ID: $(c.id)) and state (ID: $(s.id))")
        new{T}(c, s)
    end
end

export Fragment

"""
    Fragment 

A [`Fragment`](@ref) is a type overload for `Pose{Segment}` and therefore does not
contain a root/origin. These are usually used as temporary carriers of
information, without the ability to be directly incorporated in simulations.


# See also
[`fragment`](@ref)
"""
const Fragment = Pose{Segment}

Pose(::Type{T}, frag::Fragment) where {T <: AbstractFloat} = begin
    top = Topology(frag.graph.name, 1)
    state = State{T}()
    state.id = top.id
    pose = Pose(top, state)

    # Note: Frag gets copied in append_fragment_as_new_segment
    ProtoSyn.append_fragment_as_new_segment!(pose, frag)

    ProtoSyn.request_i2c!(state; all=true)
    return pose
end

Pose(frag::Fragment) = Pose(Units.defaultFloat, frag)

Pose() = Pose(Topology("unk", -1), State())

function reindex(pose::Pose)
    reindex(pose.graph)
    reindex(pose.state)
end