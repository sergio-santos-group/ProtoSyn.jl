"""
    append!(pose::Pose{Topology}, frag::Fragment)

Append a [Fragment](@ref) `frag` as a new [`Segment`](@ref) to the given
[Pose](@ref) `pose`. This function overwrites `pose`.

*This function is a Base module overload.*

# Examples
```jldoctest
julia> append!(pose, frag)
Pose{Topology}(Topology{/2a3d:35776}, State{Float64}:
 Size: 2280
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
Base.append!(pose::Pose{Topology}, frag::Fragment) = begin
    # This function is called by the `Builder.build` function.

    !isfragment(frag) && error("Invalid fragment")
    _frag = copy(frag)
    
    # Merge the fragment graph (Segment) to the pose graph (Topology).
    push!(pose.graph, _frag.graph)

    # Merge the fragment state to the pose state.
    Base.append!(pose.state, _frag.state)
    
    # Make sure the fragment graph has the same origin of the new pose.
    root_residue = origin(_frag.graph).container
    setparent!(origin(_frag.graph), root(pose.graph))
    setparent!(root_residue, root(pose.graph).container)

    # Re-index the pose to account for the new segment/residue/atoms
    reindex(pose.graph)
    pose
end

