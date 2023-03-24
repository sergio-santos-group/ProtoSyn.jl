export fragment
"""
    fragment(pose::Pose{Topology})
    
Return a [`Fragment`](@ref) from a given [Pose](@ref pose-types) `pose`. The pose must have a
single [`Segment`](@ref).

    fragment(pose::Pose{Topology}, selection::ProtoSyn.AbstractSelection)

Return a [`Fragment`](@ref) from a list of residues retrieved from the given
`selection` when applied to the provided [Pose](@ref pose-types) `pose`. If not yet of
selection type [`Residue`](@ref), the `selection` will be promoted to
[`Residue`](@ref) selection type (with the default `any` aggregating function).
The resulting list of residues must be contiguous (a connected graph of
[`Residue`](@ref) instances parenthoods). These will constitute the unique
[`Segment`](@ref) of the resulting [`Fragment`](@ref).

!!! ukw "Note:"
    A [`Fragment`](@ref) is a `Pose{Segment}`, without a root/origin. These are
    usually used as temporary carriers of information, without the ability to be
    directly incorporated in simulations.

# Examples
```
julia> frag = fragment(pose)
Fragment(Segment{/UNK:9547}, State{Float64}:
 Size: 343
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)

julia> frag = fragment(pose, rid"1:10")
Fragment(Segment{/UNK:58266}, State{Float64}:
 Size: 160
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function fragment(pose::Pose{Topology})
    length(pose.graph) != 1 && error("Only topologies with a single Segment can be turned into fragments")
    
    topology = copy(pose.graph)
    segment = topology[1]
    state = splice!(copy(pose.state), 1:count_atoms(segment))
    
    # Detach segment from the old root. This includes removing any parenthood to
    # the origin on any atom or residue.
    detach(segment)

    segment.id = state.id = genid()
    segment.name = topology.name
    segment.container = nothing

    Pose(segment, state)
end


function fragment(pose::Pose{Topology}, selection::ProtoSyn.AbstractSelection)
    # Assumes all residues selected belong to the same Segment

    sele = promote(selection, Residue)
    if count(sele(pose)) === 0
        @error "The provided selection yielded no Residue instances for fragmentation."
        return nothing
    end
    if !ProtoSyn.is_contiguous(pose, sele)
        error("Tried to fragment a non-contigous selection of residues.")
    end

    @assert length(unique([res.container.id for res in sele(pose, gather = true)])) == 1 "Tried to fragment a contiguous selection of residues belonging to different Segments."

    # Get a copy of the selected residues as a new Segment
    residues       = sele(pose, gather = true)
    copied_segment = copy(residues[1].container)
    residues       = sele(copied_segment, gather = true)
    segment        = Segment(residues[1].container.name, 1)
    segment.items  = residues
    segment.size   = length(segment.items)
    segment.code   = segment.name[1]

    # Pop parent of any Atom instance connected outside the new fragment (and
    # corresponding Residue container)
    atoms = collect(eachatom(segment))
    for atom in atoms
        for bond in atom.bonds
            !(bond in atoms) && begin

                # Unbond
                i = findfirst(atom, bond.bonds)
                i !== nothing && deleteat!(bond.bonds, i)
                
                j = findfirst(bond, atom.bonds)
                j !== nothing && deleteat!(atom.bonds, j)

                # Remove parenthood (Atom level)
                ProtoSyn.isparent(bond, atom) && begin
                    ProtoSyn.popparent!(atom)
                end
                ProtoSyn.isparent(atom, bond) && begin
                    ProtoSyn.popparent!(bond)
                end

                # Remove parenthood (Residue level)
                ProtoSyn.isparent(bond.container, atom.container) && begin
                    ProtoSyn.popparent!(atom.container) 
                end
                ProtoSyn.isparent(atom.container, bond.container) && begin
                    ProtoSyn.popparent!(bond.container) 
                end
            end
        end
    end

    # Set the residue container
    for residue in segment.items
        residue.container = segment
    end

    # Get a copy of the selected residues' State
    n_atoms = count_atoms(segment)
    state = ProtoSyn.State(n_atoms)
    for (index, atom) in enumerate(eachatom(segment))
        state.items[index + state.index_offset]        = copy(pose.state[atom])
        state.items[index + state.index_offset].parent = state
        state.items[index + state.index_offset].index  = index
        state.x[:, index] = copy(pose.state.x[:, atom.index])
    end

    segment.id = state.id = genid()
    reindex(segment)
    reindex(state)

    return Pose(segment, state)
end

# function fragment(coords::Vector{Vector{T}}) where {T <: AbstractFloat}
#     N = length(coords)
#     state = State(N)
#     segment = Segment("UNK", -1)
#     for i in 1:N
#         res = Residue!(segment, "UNK", i)
#         state[i].t = coords[i]
#         Atom!(res, "X$i", i, i, "X")
#     end

#     return Pose(segment, state)
# end


"""
    fragment!(pose::Pose{Topology}, selection::ProtoSyn.AbstractSelection; [keep_downstream_position::Bool = true])

Return a [`Fragment`](@ref) from a list of residues retrieved from the given
`AbstractSelection` `selection` when applied to the provided [Pose](@ref pose-types)
`pose`. If not yet of selection type [`Residue`](@ref), the `selection` will be
promoted to [`Residue`](@ref) selection type (with the default `any` aggregating
function). The resulting list of residues must be contiguous (a connected graph
of [`Residue`](@ref) instances parenthoods). These will constitute the unique
[`Segment`](@ref) of the resulting [`Fragment`](@ref). **In opposition to the
[`fragment`](@ref) method, this function will remove the fragmented
[`Residue`](@ref) instances from the original [`Pose`](@ref) (using the
[`pop_residue!`](@ref) method).** If `keep_downstream_position` is set to `true`
(is, by default), the downstream [`Residue`](@ref) position is maintained (by
calling [`request_c2i!`](@ref) and [`sync!`](@ref) methods).

!!! ukw "Note:"
    A [`Fragment`](@ref) is a `Pose{Segment}`, without a root/origin. These are
    usually used as temporary carriers of information, without the ability to be
    directly incorporated in simulations.

# Examples
```
julia> frag = fragment!(pose, rid"19:26")
Fragment(Segment{/UNK:9547}, State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function fragment!(pose::Pose{Topology}, selection::ProtoSyn.AbstractSelection;
    keep_downstream_position::Bool = true)
    
    frag = ProtoSyn.fragment(pose, selection)
    residues = promote(selection, Residue)(pose, gather = true)
    for residue in reverse(residues)
        ProtoSyn.pop_residue!(pose, residue,
            keep_downstream_position = keep_downstream_position)
    end

    return frag
end


export isfragment
"""
    isfragment(pose::Pose)

Return `true` if the given `pose` [Graph](@ref state-types) is a single non-empty
[`Segment`](@ref) (with no container).

# See also
[`fragment`](@ref)

# Examples
```jldoctest
julia> isfragment(frag)
true

julia> isfragment(pose)
false
```
"""
@inline isfragment(pose::Pose) = begin
    return !(hascontainer(pose.graph) || isempty(pose.graph)) && typeof(pose.graph) === Segment
end


"""
    append_fragment_as_new_segment!(pose::Pose{Topology}, frag::Fragment)

Append a [`Fragment`](@ref) `frag` as a new [`Segment`](@ref) to the given
[`Pose`](@ref) `pose`. This function overwrites `pose`.

# See also
[`isfragment`](@ref) [`fragment`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.append_fragment_as_new_segment!(pose, frag)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 373
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function append_fragment_as_new_segment!(pose::Pose{Topology}, frag::Fragment)
    # This function is called by the `ProtoSyn.build` function.

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
    reindex(pose.graph, set_ascendents = true)
    pose
end


export append_fragment!
"""
    append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, frag::Pose{Segment}; op = "α")

Add the [`Fragment`](@ref) `frag` to the given [`Pose`](@ref) `pose`, appending
it after the given [`Residue`](@ref) `residue`. This residue and the new
[`Fragment`](@ref) `frag` will be connected using operation `op` ("α" by
default) of the given [`LGrammar`](@ref) `grammar`. Request internal to
cartesian coordinate conversion and return the altered [`Pose`](@ref) `pose`.

# See also
[`insert_fragment!`](@ref ProtoSyn.insert_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Pose{Segment}; ::Any))

# Examples
```jldoctest
julia> ProtoSyn.append_fragment!(pose, pose.graph[1][end], res_lib, frag)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 373
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, frag::Pose{Segment}; op = "α")

    sync!(pose)

    # Soft `uncap!` if is C terminal. Might change in future versions.
    if ProtoSyn.Peptides.is_C_terminal(residue) && residue["OXT"] !== nothing
        ProtoSyn.pop_atom!(pose, residue["OXT"])
    end

    # Insert the fragment residues in the pose.graph and set
    # frag_residue.container (of each residue in the fragment) to be the segment
    # of the "parent" residue (automatically on `insert!`)
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
    
    ProtoSyn.request_i2c!(pose.state; all = true)
    return pose
end


export insert_fragment!
"""
    insert_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, frag::Pose{Segment}; op = "α")

Insert the [`Fragment`](@ref) `frag` in the given `pose`, on the position of the
provided [`Residue`](@ref) instance `residue` (the `residue` gets shifted
downstream). This first downstream [`Residue`](@ref) and the new
[`Fragment`](@ref) will be connected using operation `op` ("α" by default) from
[`LGrammar`] `grammar`. Also connects to the upstream [`Residue`](@ref)
instance, using the same operation. Request internal to cartesian coordinate
conversion and return the altered [`Pose`](@ref) `pose`.

!!! ukw "Note:"
    Consider using more specific versions of this function, see [`Peptides.insert_fragment!`](@ref Peptides.insert_fragment!)

# See also
[`append_fragment!`](@ref ProtoSyn.append_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Pose{Segment}; ::Any))

# Examples
```jldoctest
julia> ProtoSyn.insert_fragment!(pose, pose.graph[1][1], res_lib, frag)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 373
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function insert_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, frag::Pose{Segment}; op = "α")

    # Setup
    sync!(pose)
    r_index         = residue.index
    frag_size       = length(eachresidue(frag.graph))
    anchor_id       = residue.parent.id
    connect_to_root = false
    if residue.parent == ProtoSyn.root(pose.graph).container
        connect_to_root = true
    end

    # Remove any bonds to the parent residue (if not connected to root)
    !connect_to_root && begin
        for atom in residue.items
            for bond in atom.bonds
                if bond in residue.parent.items
                    @info "Unbonding upstream connection ... $atom - $bond"
                    ProtoSyn.unbond!(pose, atom, bond,
                        keep_downstream_position = false)
                end
            end
        end
    end

    # Insert the fragment residues in the pose.graph and set
    # frag_residue.container (of each residue in the fragment) to be the segment
    # of the "parent" residue (automatically on `insert!`)
    @info "Shifting $residue downstream by $frag_size residues ..."
    insert!(residue.container, residue.index, frag.graph.items)
    insert!(pose.state, residue.items[1].index, frag.state)

    # Perform link operation. Requires correct indexes. Set ascendents is set to
    # false because some residues are still orphan (would cause error/bug). Set:
    # - Distance/angle/dihedrals in the fragment first residue
    # - Parent/children in the newly bonded atoms
    # - Parent/children in the newly bonded residues
    # - Atom bonds in the newly bonded atoms
    reindex(pose.graph, set_ascendents = false)
    reindex(pose.state)

    if connect_to_root
        # Case needs to be connected to root
        _root = ProtoSyn.root(pose.graph)
        @info "Joining upstream $(_root.container) - $(pose.graph[1][r_index]) ..."
        for atom in pose.graph[1][r_index].items
            atom.parent === nothing && ProtoSyn.setparent!(atom, _root)
        end
        ProtoSyn.setparent!(pose.graph[1][r_index], _root.container)
    else
        anchor = residue.container[anchor_id]
        @info "Joining upstream $anchor - $(pose.graph[1][r_index]) ..."
        grammar.operators[op](anchor, pose, residue_index = r_index)
    end

    @info "Connecting downstream ..."
    anchor = residue.container[r_index + frag_size - 1]
    
    # Remove all inter-residue parenthood relationships to upstream residue
    for atom in residue.items
        if !(atom.parent in residue.items)
            @info "Removing downstream parenthoods ... $atom & $residue"
            ProtoSyn.popparent!(atom)
            ProtoSyn.popparent!(residue)
        end
    end

    @info "Joining downstream $anchor - $(pose.graph[1][residue.index])"
    grammar.operators[op](anchor, pose, residue_index = residue.index)

    # Set correct ascedents
    reindex(pose.graph)

    ProtoSyn.request_i2c!(pose.state, all = true)
    return pose
end