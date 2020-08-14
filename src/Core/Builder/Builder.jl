module Builder

using YAML

using ..ProtoSyn
using ..ProtoSyn.Units: tonumber

#region fragment ----------------------------------------------------------------

isfragment(p::Pose) = !(hascontainer(p.graph) || isempty(p.graph))


"""
    fragment([T=Float64], grammar::LGrammar, derivation) where {T <: AbstractFloat}
    
Return a new fragment (`Pose` instance with just a single `Segment`) using the
given `derivation` sequence on the provided `grammar` instructions. Note: A 
fragment does not contain a `Topology` instance.

# Examples
```jldoctest
julia> frag = fragment(Float64, reslib, seq"AAA")

julia> frag = fragment(reslib, seq"AAA")
```
"""
function fragment(::Type{T}, grammar::LGrammar, derivation) where {T <: AbstractFloat}

    state = State{T}()
    seg = Segment("UNK", 1)
    seg.code = 'A'

    stack = Fragment[]
    opstack = Function[]
    parent::Opt{Fragment} = nothing

    for letter in derivation
        if isop(grammar, letter)
            op = getop(grammar, letter)
            push!(opstack, op)
        elseif letter == '['
            push!(stack, parent)
        elseif letter == ']'
            parent = pop!(stack)
        elseif isvar(grammar, letter)
            frag = getvar(grammar, letter)
            
            frag2 = copy(frag)

            push!(seg, frag2.graph.items...) # Appending the residues to the segment
            append!(state, frag2.state)      # Merging the 2 states
            if parent !== nothing
                join = isempty(opstack) ? grammar.defop : pop!(opstack)
                join(parent.graph[end], frag2) # Adding ascendents and bonds correctly
            end
            parent = frag2
        end
    end
    reindex(seg)
    seg.id = state.id = ProtoSyn.genid()

    return Pose(seg, state)
end

fragment(grammar::LGrammar, derivation) = fragment(Units.defaultFloat, grammar, derivation)


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
    #(imin,imax) = extrema(map(at->at.index, eachatom(segment)))
    #state = splice!(pose.state, imin:imax)
    state = splice!(pose.state, 1:count_atoms(segment))
    detach(segment)
    segment.id = state.id = genid()
    segment.name = topology.name

    Pose(segment, state)
end

Base.detach(s::Segment) = begin
    root = origin(s)
    for at in root.children
        isparent(root, at) && popparent!(at)
    end
    hascontainer(s) && delete!(s.container, s)
end

#endregion fragment


"""
    build([T=Float64,] grammar::LGrammar, derivation)

Build a `Pose{Topology}` using the given `derivation` sequence on the provided
`grammar` instructions.

"""
function build(::Type{T}, grammar::LGrammar, derivation) where {T<:AbstractFloat}
    top = Topology("UNK", 1)
    state = State{T}()
    state.id = top.id
    pose = Pose(top, state)

    if !isempty(derivation)
        frag = fragment(T, grammar, derivation)
        append!(pose, frag) # Appending the fragment (which is a segment) to the Topology
        
        ProtoSyn.request_i2c(state; all=true)
    end
    pose
end
build(grammar::LGrammar, derivation) = build(Float64, grammar, derivation)


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


"""
    append!(pose::Pose{Topology}, frag::Fragment)

Append a fragment as a new segment.
Note: This function is called by the `build` function.
"""
Base.append!(pose::Pose{Topology}, frag::Fragment) = begin

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
the given `pose`, inserting it ON tHe POSITION of the given `residue`. This
residue and the new fragment will be connected using operation `op` ("α" by
default). only downstream of the insertion. If `connect_upstream` is set to true
(is, by default), also connect to the upstream residues. Return the altered
`pose`.

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
    insert!(residue.container, residue.index, frag.graph.items)

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


"""
    pop!(pose::Pose{Topology}, atom::Atom)::Pose{Atom}

Pop and return the desired `atom` from the given `pose`. In order to do this,
perform the following actions:
 - Unset parents/children
 - Unbond neighbours
 - Remove from graph
 - Remove from state
 - Set new ascendents
 - Update the container itemsbyname

# Examples
```jldoctest
julia> pop!(pose, pose.graph[1][1][2])
```
"""
function Base.pop!(pose::Pose{Topology}, atom::Atom)::Pose{Atom}

    if atom.container.container.container !== pose.graph
        error("Given Atom does not belong to the provided topology.")
    end

    # Save information to return 
    popped_atom = Atom(atom.name, 1, 1, atom.symbol)

    # Save children and parent of this atom (will be removed in next step)
    children    = copy(atom.children)
    grandparent = atom.parent

    # Unset parents/children and unbond neighbours
    for i = length(atom.bonds):-1:1   # Note the reverse loop
        other = atom.bonds[i]
        unbond(pose, atom, other)
    end

    # Using saved children and parent, set all child.parent to be this
    # atom.parent
    for child in children
        child.parent = grandparent
    end

    # Remove from graph
    deleteat!(atom.container.items, findfirst(atom, atom.container.items))
    atom.container.size -= 1

    # Remove from state
    popped_state = splice!(pose.state, atom.index)

    # Reindex and set ascendents
    reindex(pose.graph)

    # Update container 'itemsbyname'
    pop!(atom.container.itemsbyname, atom.name)

    # Set common ID
    popped_atom.id = popped_state.id = genid()

    return Pose(popped_atom, popped_state)
end


export setdihedral!

"""
    setdihedral!(s::State, at::Atom, val::T) where {T <: AbstractFloat}

Set the dihedral in `Atom` `at` of `State` `s` to be `val` (in radians)

# Examples
```jldoctest
julia> setdihedral!(pose.state, pose.graph[1][1][end], π)
```
"""
@inline setdihedral!(s::State, at::Atom, val::T) where {T <: AbstractFloat} = begin
    s[at].Δϕ = val - s[at].ϕ
    ProtoSyn.request_i2c(s)
    s
end

end