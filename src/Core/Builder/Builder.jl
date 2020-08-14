module Builder

using YAML

using ..ProtoSyn
using ..ProtoSyn.Units: tonumber

#region fragment ----------------------------------------------------------------

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


function append_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op = "α")

    frag = Builder.fragment(grammar, derivation)

    # Inserts the fragment residues in the Pose graph and sets residue.container
    # println("Insert residues in position $(residue.index + 1)/$(length(residue.container.items))")
    insert!(residue.container, residue.index + 1, frag.graph.items)

    # Inserts the fragment atoms state in the pose.state
    # println("Insert residues in position $(residue.items[end].index + 1)/$(length(pose.state.items)-3)")
    insert!(pose.state, residue.items[end].index + 1, frag.state)

    # Sets:
    # - Distance/angle/dihedrals in the fragment first residue
    # - Parent/children in the newly bonded atoms
    # - Parent/children in the newly bonded residues
    # - Atom bonds in the newly bonded atoms

    reindex(pose.graph, set_ascendents = false)

    grammar.operators[op](residue, pose, residue_index = residue.index + 1)

    # Reindex and define new ascendents
    reindex(pose.graph)
    ProtoSyn.request_i2c(pose.state; all=true)
end


function insert_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op = "α", connect_upstream = true)

    residue_index = residue.index
    frag = Builder.fragment(grammar, derivation)

    insert!(residue.container, residue.index, frag.graph.items)
    insert!(pose.state, residue.items[1].index, frag.state)

    # Remove all references of the origin as the parent of any of its children
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

    # Reindex to account for the inserted residues
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
end


# function mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op = "α")

#     @assert length(derivation) == 1 "Derivation must have length = 1."
#     residue_index = residue.index

#     insert_residues!(pose, residue, grammar, derivation; op = op)
#     # pop!(pose, residue.container[residue.index])



#     # bond(pose, residue.container[residue_index], residue.container[residue_index + 1], grammar, op = op)
# end

end