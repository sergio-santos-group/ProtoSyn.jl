"""
    append!(pose::Pose{Topology}, frag::Fragment)

Append a fragment as a new segment.
Note: This function is called by the `build` function.
"""
Base.append!(pose::Pose{Topology}, frag::Fragment) = begin

    println(frag, " ", isfragment(frag))
    !isfragment(frag) && error("invalid fragment")
    
    # Merge the fragment graph (Segment) to the pose graph (Topology).
    push!(pose.graph, frag.graph)

    # Merge the fragment state to the pose state.
    Base.append!(pose.state, frag.state)
    
    # Make sure the fragment graph has the same origin of the new pose.
    root_residue = root(frag.graph).container
    setparent!(root(frag.graph), origin(pose.graph))

    setparent!(root_residue, origin(pose.graph).container)

    # Re-index the pose to account for the new segment/residue/atoms
    reindex(pose.graph)
    pose
end

