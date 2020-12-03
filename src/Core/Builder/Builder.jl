module Builder

using YAML

using ..ProtoSyn
using ..ProtoSyn.Units: tonumber

include("base.jl")

"""
    build(grammar::LGrammar{T}, derivation)

Build a `Pose{Topology}` using the given `derivation` sequence on the provided
`grammar` instructions.

"""
function build(grammar::LGrammar{T}, derivation) where {T<:AbstractFloat}
    top = Topology("UNK", 1)
    state = State{T}()
    state.id = top.id
    pose = Pose(top, state)

    if !isempty(derivation)
        frag = fragment(grammar, derivation)
        append!(pose, frag) # Appending the fragment (which is a segment) to the Topology
        
        ProtoSyn.request_i2c(state; all=true)
    end
    pose
end
# build(grammar::LGrammar, derivation) = build(Float64, grammar, derivation)


export @seq_str

"""
    @seq_str -> Vector{String}

Construct a vector of strings from the provided string. 

# Examples
```jldoctest
julia> seq"ABC"
3-element Array{String,1}:
 "A"
 "B"
 "C"
```
"""
macro seq_str(s); [string(c) for c in s]; end


# """
#     append!(pose::Pose{Topology}, frag::Fragment)

# Append a fragment as a new segment.
# Note: This function is called by the `build` function.
# """
# Base.append!(pose::Pose{Topology}, frag::Fragment) = begin

#     !isfragment(frag) && error("invalid fragment")
    
#     # Merge the fragment graph (Segment) to the pose graph (Topology).
#     push!(pose.graph, frag.graph)

#     # Merge the fragment state to the pose state.
#     Base.append!(pose.state, frag.state)
    
#     # Make sure the fragment graph has the same origin of the new pose.
#     root_residue = root(frag.graph).container
#     setparent!(root(frag.graph), origin(pose.graph))

#     setparent!(root_residue, origin(pose.graph).container)

#     # Re-index the pose to account for the new segment/residue/atoms
#     reindex(pose.graph)
#     pose
# end


"""
    append_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op = "α")

Based on the provided `grammar`, add the residue sequence from `derivation` to
the given `pose`, appending it AFTER the given `residue`. This residue and the
new fragment will be connected using operation `op` ("α" by default). Return the
altered `pose`.

# Examples
```jldoctest
julia> append_residues!(pose, pose.graph[1][1], reslib, seq"A")
```
"""
function append_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op = "α")

    # Build the fragment to append
    frag = Builder.fragment(grammar, derivation)

    # Insert the fragment residues in the pose.graph and set
    # frag_residue.container (of each residue in the fragment) to be the segment
    # of the "parent" residue
    insert!(residue.container, residue.index + 1, frag.graph.items)

    # Inserts the fragment atoms state in the pose.state
    insert!(pose.state, residue.items[end].index + 1, frag.state)

    # Perform link operation. Requires correct indexes. Set ascendents is set to
    # false because some residues are still orphan (would cause error/bug). Set:
    # - Distance/angle/dihedrals in the fragment first residue
    # - Parent/children in the newly bonded atoms
    # - Parent/children in the newly bonded residues
    # - Atom bonds in the newly bonded atoms
    reindex(pose.graph, set_ascendents = false)
    grammar.operators[op](residue, pose, residue_index = residue.index + 1)

    # Reindex to define new ascendents
    reindex(pose.graph)
    
    ProtoSyn.request_i2c(pose.state; all=true)
    return pose
end


"""
    insert_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op = "α", connect_upstream = true)

Based on the provided `grammar`, add the residue sequence from `derivation` to
the given `pose`, inserting it ON THE POSITION of the given `residue` (the
`residue` gets shifted downstream). The first downstream residue and the new
fragment will be connected using operation `op` ("α" by default). If
`connect_upstream` is set to true (is, by default), also connect to the upstream
residues. Return the altered `pose`.

# Examples
```jldoctest
julia> insert_residues!(pose, pose.graph[1][2], reslib, seq"A")
```
"""
function insert_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op = "α", connect_upstream = true)

    residue_index = residue.index
    frag = Builder.fragment(grammar, derivation)

    # Insert the fragment residues in the pose.graph and set
    # frag_residue.container (of each residue in the fragment) to be the segment
    # of the "parent" residue
    ProtoSyn.insert!(residue.container, residue.index, frag.graph.items)

    # Inserts the fragment atoms state in the pose.state
    insert!(pose.state, residue.items[1].index, frag.state)

    # Remove all references of the origin as the parent of any of its children
    # (if said children belong to the residue that is being displaced)
    # and remove said children from origin.children
    pose_origin = ProtoSyn.origin(pose.graph) # This is an Atom
    parent_is_origin = residue.parent == pose_origin.container
    if parent_is_origin
        for child in pose_origin.children
            child in residue.items && popparent!(child)
        end
        
        # Add the first atom of the appendage as a child of origin AND set its
        # parent as being the origin
        _root = ProtoSyn.root(frag.graph)
        setparent!(_root, pose_origin)
        setparent!(_root.container, pose_origin.container)
        popparent!(residue)
    end

    # Perform link operation. Requires correct indexes. Set ascendents is set to
    # false because some residues are still orphan (would cause error/bug). Set:
    # - Distance/angle/dihedrals in the fragment first residue
    # - Parent/children in the newly bonded atoms
    # - Parent/children in the newly bonded residues
    # - Atom bonds in the newly bonded atoms
    reindex(pose.graph, set_ascendents = false)
    grammar.operators[op](frag.graph[end], pose, residue_index = residue_index + length(frag.graph))


    # Case we are inserting between two pre-existing residues (so far, the same
    # operation will be used in both cases. Might be useful to differentiate
    # between left and right operation.)
    if connect_upstream
        for child in pose_origin.children
            child in residue.container[residue_index].items && popparent!(child)
        end
        popparent!(residue.container[residue_index])

        grammar.operators[op](residue.container[residue_index - 1], pose, residue_index = residue_index)
    end
    
    # Reindex to set correct ascendents
    reindex(pose.graph)
    
    ProtoSyn.request_i2c(pose.state; all=true)
    return pose
end


# """
#     pop!(pose::Pose{Topology}, atom::Atom)::Pose{Atom}

# Pop and return the desired `atom` from the given `pose`. In order to do this,
# perform the following actions:
#  - Unset parents/children
#  - Unbond neighbours
#  - Remove from graph
#  - Remove from state
#  - Set new ascendents
#  - Update the container itemsbyname

# # Examples
# ```jldoctest
# julia> pop!(pose, pose.graph[1][1][2])
# ```
# """
# function Base.pop!(pose::Pose{Topology}, atom::Atom)::Pose{Atom}

#     if atom.container.container.container !== pose.graph
#         error("Given Atom does not belong to the provided topology.")
#     end

#     # Save information to return 
#     popped_atom = Atom(atom.name, 1, 1, atom.symbol)

#     # Save children and parent of this atom (will be removed in next step)
#     children    = copy(atom.children)
#     grandparent = atom.parent

#     # Unset parents/children and unbond neighbours
#     for i = length(atom.bonds):-1:1   # Note the reverse loop
#         other = atom.bonds[i]
#         ProtoSyn.unbond(pose, atom, other)
#     end

#     # Using saved children and parent, set all child.parent to be this
#     # atom.parent
#     for child in children
#         # child.parent = grandparent
#         child.parent = origin(pose.graph)
#     end

#     # Remove from graph
#     deleteat!(atom.container.items, findfirst(atom, atom.container.items))
#     atom.container.size -= 1

#     # Remove from state
#     popped_state = splice!(pose.state, atom.index)

#     # Reindex and set ascendents
#     reindex(pose.graph)

#     # Update container 'itemsbyname'
#     pop!(atom.container.itemsbyname, atom.name)

#     # Set common ID
#     popped_atom.id = popped_state.id = genid()

#     return Pose(popped_atom, popped_state)
# end

end