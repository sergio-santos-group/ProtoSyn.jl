function Base.pop!(pose::Pose{Topology}, residue::Residue)

    for atom in residue.items
        pop!(pose, atom)
    end
end


# Base.pop!(pose::Pose{Topology}, r::Residue) = begin

#     println("Deleting residue $r")

#     if r.container.container !== pose.graph
#         error("given residue does not belong to the provided topology")
#     end

#     # detach residue from parents (bonds, parents/children)
#     # hascontainer(r) && delete!(r.container, r)
#     for atom in eachatom(r)
#         for i = length(atom.bonds):-1:1   # Note the reverse loop
#             other = atom.bonds[i]

#             in(other, r) && continue
#             Builder.unbond(pose, atom, other)
#         end
#     end

#     # remove node states from parent state and create
#     # a new state for this residue
#     deleteat!(r.container.items, findfirst(r, r.container.items))
#     r.container.size -= 1
#     for child in origin(pose.graph).children
#         child in r.items && popparent!(child) 
#     end

#     st = splice!(pose.state, r[1].index:r[end].index)
#     reindex(pose.graph) # Also sets ascendents
    
#     # new common ID
#     r.id = st.id = genid()

#     #return (res, st)
#     Pose(r, st)
# end

Base.pop!(pose::Pose{Topology}, seg::Segment) = begin

    if seg.container !== pose.graph
        error("given residue does not belong to the provided topology")
    end

    # remove node states from parent state and create
    # a new state for this segment
    st = splice!(pose.state, seg[1][1].index:seg[end][end].index)
    
    pop!(pose.graph.items, seg)
    # deleteat!(pose.graph.items, findall(x -> x == seg, pose.graph.items))
    # pose.graph.size -= 1

    reindex(pose.graph)
    
    # new common ID
    seg.id = st.id = genid()
    Pose(seg, st)
end


Base.pop!(top::Topology, seg::Segment) = begin
    deleteat!(top.items, findall(x -> x == seg, top.items))
    top.size -= 1
end