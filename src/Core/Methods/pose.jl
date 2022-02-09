function recoverfrom!(pose::Pose, backup::Pose)
    pose.state = copy(backup.state)
    pose.graph = copy(backup.graph)
end


"""
    merge(pose1::Pose, pose2::Pose)

Merge the two given poses, creating a new `Pose` in the process.

# Examples
```
julia> ProtoSyn.merge(pose, pose_mod)
Pose{Topology}(Topology{/merged:32083}, State{Float64}:
 Size: 686
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
...
```
"""
function merge(pose1::Pose, pose2::Pose)::Pose

    function merge_segments(pose::Pose)
        for segment in pose.graph.items
            s = copy(segment)
            push!(graph, s)
            # Set new parenthood of the first residue in the segment
            ProtoSyn.setparent!(s[1], root.container)
            # Set new parenthood of the first atom in the segment
            ProtoSyn.setparent!(s[1][1], root)
        end
    end

    # Merge graphs
    graph = Topology("merged", -1)
    root = ProtoSyn.root(graph)
    merge_segments(pose1)
    merge_segments(pose2)

    # Merge states
    state = State(pose1.state.size + pose2.state.size)
    state.x[:, 1:pose1.state.size] = pose1.state.x.coords
    state.x[:, (pose1.state.size+1):end] = pose2.state.x.coords

    reindex(graph, set_ascendents = true)
    reindex(state)
    graph.id = state.id = genid()
    sync!(state, graph)
    return Pose(graph, state)
end


"""
    merge!(pose1::Pose, pose2::Pose)

Merge the two given poses, updating/overwritting the given `pose1`.

# Examples
```
julia> ProtoSyn.merge!(pose, pose_mod)
Pose{Topology}(Topology{/merged:10313}, State{Float64}:
 Size: 748
 i2c: false | c2i: true
 Energy: Dict(:Total => Inf)
)
...
```
"""
function merge!(pose1::Pose, pose2::Pose)

    # Merge graphs
    root = ProtoSyn.root(pose1.graph)
    for segment in pose2.graph.items
        s = copy(segment)
        push!(pose1.graph, s)
        # Set new parenthood of the first residue in the segment
        ProtoSyn.setparent!(s[1], root.container)
        # Set new parenthood of the first atom in the segment
        ProtoSyn.setparent!(s[1][1], root)
    end

    # Merge states (including forces)
    pose1.state.f = hcat(pose1.state.f, pose2.state.f)
    pose1.state.x.coords = hcat(pose1.state.x.coords, pose2.state.x.coords)
    for item in pose2.state.items[4:end]
        item.parent = pose1.state
    end
    pose1.state.items = vcat(pose1.state.items, pose2.state.items[4:end])
    pose1.state.size = length(pose1.state.items) - 3

    reindex(pose1.graph, set_ascendents = true)
    reindex(pose1.state)
    return pose1
end


"""
    symexp(pose::Pose, reps::Vector{Int}, unit_cell_dims::Vector{T}) where {T <: AbstractFloat}

Return a symmetry expanded [Pose](@ref). Create N copies of the given `pose` in
all 3 symmetry axis of a cubic lattice, where `reps` is the number of copies in
each of the dimensions X, Y and Z (N is, therefore, reps[1]*reps[2]*reps[3]).
Length of `reps` must be 3. `unit_cell_dims` sets the distance in each of
dimension to translate the copies, in Angstrom Å. Length of `unit_cell_dims`
must be 3. Copies the given `pose`, returning a new struct.

# See also
`symexp!` `merge`

# Examples
```
julia> ProtoSyn.symexp(pose, [2, 2, 2], [50.0, 50.0, 50.0])
Pose{Topology}(Topology{/UNK:59312}, State{Float64}:
 Size: 9261
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
...
```
"""
function symexp(pose::Pose, reps::Vector{Int}, unit_cell_dims::Vector{T}) where {T <: AbstractFloat}
    @assert length(reps)==3 "`reps` should be a Vector{Int} with 3 numbers - x, y and z number of repetitions"
    @assert length(unit_cell_dims)==3 "`unit_cell_dims` should be a Vector{AbstractFloat} with 3 numbers - x, y and z distances of the bounding box"

    _pose = copy(pose)
    return symexp!(_pose, reps, unit_cell_dims)
end


"""
    symexp!(pose::Pose, reps::Vector{Int}, unit_cell_dims::Vector{T}) where {T <: AbstractFloat}

Return a symmetry expanded [Pose](@ref). Create N copies of the given `pose` in
all 3 symmetry axis of a cubic lattice, where `reps` is the number of copies in
each of the dimensions X, Y and Z (N is, therefore, reps[1]*reps[2]*reps[3]).
Length of `reps` must be 3. `unit_cell_dims` sets the distance in each of
dimension to translate the copies, in Angstrom Å. Length of `unit_cell_dims`
must be 3. Copies the given `pose`, returning a new struct. Updates/overwrites
the given `pose`.

# See also
`symexp` `merge`

# Examples
```
julia> ProtoSyn.symexp!(pose, [2, 2, 2], [50.0, 50.0, 50.0])
Pose{Topology}(Topology{/merged:10313}, State{Float64}:
 Size: 748
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
...
```
"""
function symexp!(pose::Pose, reps::Vector{Int}, unit_cell_dims::Vector{T}) where {T <: AbstractFloat}
    @assert length(reps)==3 "`reps` should be a Vector{Int} with 3 numbers - x, y and z number of repetitions"
    @assert length(unit_cell_dims)==3 "`unit_cell_dims` should be a Vector{AbstractFloat} with 3 numbers - x, y and z distances of the bounding box"

    x = unit_cell_dims[1]
    y = unit_cell_dims[2]
    z = unit_cell_dims[3]

    i_pose = copy(pose) # pose without the newly added pieces
    for i in 0:reps[1]
        for j in 0:reps[2]
            for k in 0:reps[3]
                _pose = copy(i_pose)
                i==0 && j == 0 && k == 0 && continue
                translation = [i * x, j * y, k * z]
                for i in 1:_pose.state.size
                    _pose.state.x[:, i] = pose.state.x[:, i] .+ translation
                end
                ProtoSyn.request_c2i!(_pose.state)
                sync!(_pose)
                ProtoSyn.merge!(pose, _pose)    
            end # k for
        end # j for
    end # i for
    
    return pose
end # function


export fragment
"""
    fragment(pose::Pose{Topology})
    
Return a [Fragment](@ref) from a given [Pose](@ref) `pose`. The pose must have a
single [`Segment`](@ref).

    fragment(pose::Pose{Topology}, selection::ProtoSyn.AbstractSelection)

Return a [Fragment](@ref) from a list of residues retrieved from the given
`selection` when applied to the provided [Pose](@ref) `pose`. If not yet of
selection type [`Residue`](@ref), the `selection` will be promoted to
[`Residue`](@ref) selection type (with the default `any` aggregating function).
The resulting list of residues must be contiguous (a connected graph of
[`Residue`](@ref) instances parenthoods). These will constitute the unique
[`Segment`](@ref) of the resulting [Fragment](@ref).

!!! ukw "Note:"
    A [Fragment](@ref) is a `Pose{Segment}`, without a root/origin. These are
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

                # Remove parenthood
                ProtoSyn.isparent(bond.container, atom.container) && begin
                    ProtoSyn.popparent!(atom.container) 
                end
                ProtoSyn.isparent(atom.container, bond.container) && begin
                    ProtoSyn.popparent!(bond.container) 
                end
               
            end
        end
    end

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


"""
    fragment(coords::Vector{Vector{T}}) where {T <: AbstractFloat}

Return a [Fragment](@ref) from a list of `coords`. Each coordinate creates a new
[Residue](@ref) instance with a single [Atom](@ref) instance. Doesn't set
parenthoods.

# Examples
```
julia> frag = fragment(coords)
Fragment(Segment{/UNK:-1}, State{Float64}:
 Size: 343
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function fragment(coords::Vector{Vector{T}}) where {T <: AbstractFloat}
    N = length(coords)
    state = State(N)
    segment = Segment("UNK", -1)
    for i in 1:N
        res = Residue!(segment, "UNK", i)
        state[i].t = coords[i]
        Atom!(res, "X$i", i, i, "X")
    end

    return Pose(segment, state)
end


"""
    fragment!(pose::Pose{Topology}, selection::ProtoSyn.AbstractSelection; [keep_downstream_position::Bool = true])

Return a [Fragment](@ref) from a list of residues retrieved from the given
`AbstractSelection` `selection` when applied to the provided [Pose](@ref)
`pose`. If not yet of selection type [`Residue`](@ref), the `selection` will be
promoted to [`Residue`](@ref) selection type (with the default `any` aggregating
function). The resulting list of residues must be contiguous (a connected graph
of [`Residue`](@ref) instances parenthoods). These will constitute the unique
[`Segment`](@ref) of the resulting [Fragment](@ref). **In opposition to the
[`fragment`](@ref) method, this function will remove the fragmented
[`Residue`](@ref) instances from the original [`Pose`](@ref) (using the
[`pop_residue!`](@ref) method).** If `keep_downstream_position` is set to `true`
(is, by default), the downstream [`Residue`](@ref) position is maintained (by
calling [`request_c2i!`](@ref) and [`sync!`](@ref) methods).

!!! ukw "Note:"
    A [Fragment](@ref) is a `Pose{Segment}`, without a root/origin. These are
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

Append a [Fragment](@ref) `frag` as a new [`Segment`](@ref) to the given
[Pose](@ref) `pose`. This function overwrites `pose`.

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
    reindex(pose.graph)
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
                    ProtoSyn.verbose.mode && @info "Unbonding upstream connection ... $atom - $bond"
                    ProtoSyn.unbond!(pose, atom, bond,
                        keep_downstream_position = false)
                end
            end
        end
    end

    # Insert the fragment residues in the pose.graph and set
    # frag_residue.container (of each residue in the fragment) to be the segment
    # of the "parent" residue (automatically on `insert!`)
    ProtoSyn.verbose.mode && @info "Shifting $residue downstream by $frag_size residues ..."
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
        ProtoSyn.verbose.mode && @info "Joining upstream $(_root.container) - $(pose.graph[1][r_index]) ..."
        for atom in pose.graph[1][r_index].items
            atom.parent === nothing && ProtoSyn.setparent!(atom, _root)
        end
        ProtoSyn.setparent!(pose.graph[1][r_index], _root.container)
    else
        anchor = residue.container[anchor_id]
        ProtoSyn.verbose.mode && @info "Joining upstream $anchor - $(pose.graph[1][r_index]) ..."
        grammar.operators[op](anchor, pose, residue_index = r_index)
    end

    ProtoSyn.verbose.mode && @info "Connecting downstream ..."
    anchor = residue.container[r_index + frag_size - 1]
    
    # Remove all inter-residue parenthood relationships to upstream residue
    for atom in residue.items
        if !(atom.parent in residue.items)
            ProtoSyn.verbose.mode && @info "Removing downstream parenthoods ... $atom & $residue"
            ProtoSyn.popparent!(atom)
            ProtoSyn.popparent!(residue)
        end
    end

    ProtoSyn.verbose.mode && @info "Joining downstream $anchor - $(pose.graph[1][residue.index])"
    grammar.operators[op](anchor, pose, residue_index = residue.index)

    # Set correct ascedents
    reindex(pose.graph)

    ProtoSyn.request_i2c!(pose.state, all = true)
    return pose
end


"""
    insert_atom_as_children!(pose::Pose, parent_atom::Atom, atom::Atom, [atomstate::Opt{AtomState} = nothing])

Add the given [`Atom`](@ref) `atom` to the [`Pose`](@ref) `pose` graph, as a
child of `parent_atom`. Correctly sets `atom.container`, `container.size`,
`container.items_by_name`, parenthood relationships, bonds, indexes and
ascedents. If an optional [`AtomState`](@ref) `atomstate` is provided, the
inserted atom's [`State`](@ref) is set, otherwise, insert an empty
[`State`](@ref) (with all internal and cartesian coordinates set to zero).
Return the modified (in-place) [`Pose`](@ref) `pose`.

# Examples
```jldoctest
julia> ProtoSyn.insert_atom_as_children!(pose, pose.graph[1][1][1], Atom("N", 1, 1, "N"))
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 344
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function insert_atom_as_children!(pose::Pose, parent_atom::Atom, atom::Atom, atomstate::Opt{AtomState} = nothing)

    # * 1- Add atom to the pose graph
    residue = parent_atom.container
    p_atom_index_in_container = findfirst(parent_atom, residue.items)
    insert!(residue, p_atom_index_in_container + 1, atom)
    # * Note: `insert!` above already sets atom.container, container.size += 1
    # * and container.items_by_name

    # * 2- Set correct bonds and parenthood
    setparent!(atom, parent_atom)
    ProtoSyn.bond(atom, parent_atom)

    # * 3- Add atom to the pose state. Note that pose.state.items includes the
    # * origin, while pose.state.x doesn't (and therefore the index for
    # * insertion must be -3).
    p_atom_index_in_state = findfirst(pose.state[parent_atom], pose.state.items)
    if atomstate === nothing
        insert!(pose.state, p_atom_index_in_state - 2, State(1))
    else
        insert!(pose.state, p_atom_index_in_state - 2, State([atomstate]))
    end

    # * 4- Reindex graph and state, setting ascendents
    reindex(pose.state)
    reindex(pose.graph)

    return pose
end


"""
    pop_atom!(pose::Pose{Topology}, atom::Atom; [keep_downstream_position::Bool = true])

Pop and return the given [`Atom`](@ref) `atom` from the given [`Pose`](@ref)
`pose`. In order to do this, perform the following actions:
+ Unset parenthood relationships (On [`Atom`](@ref) level only);
+ Unbond neighbouring [`Atom`](@ref) instances;
+ Remove from [Graph](@ref graph-types);
+ Remove from [State](@ref state-types);
+ Set new [`ascendents`](@ref);
+ Update the `container.itemsbyname`.

If `keep_downstream_position` is set to `true` (is, by default), the
downstream [`Residue`](@ref) position is maintained (by calling
[`request_c2i!`](@ref) and [`sync!`](@ref) methods). In either case, this method
requests internal to cartesian coordinates conversion at the end (using the
[`request_i2c!`](@ref) method).

# See also
[`pop_residue!`](@ref)

# Examples
```
julia> ProtoSyn.pop_atom!(pose, pose.graph[1][1][2])
Pose{Atom}(Atom{/H:6299}, State{Float64}:
 Size: 1
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function pop_atom!(pose::Pose{Topology}, atom::Atom; keep_downstream_position::Bool = true)::Pose{Atom}

    ProtoSyn.verbose.mode && @info "Removing atom $atom ..."
    if atom.container.container.container !== pose.graph
        error("Atom $atom does not belong to the provided topology.")
    end

    # Save information to return
    popped_atom = Atom(atom.name, 1, 1, atom.symbol)

    # Unset parents/children and unbond neighbours
    # (includes children and parent)
    for i = length(atom.bonds):-1:1   # Note the reverse loop
        other = atom.bonds[i]
        # ProtoSyn.unbond already pops parenthood relationship between the two
        # atoms and sets parent of downstream atom to origin
        ProtoSyn.unbond!(pose, atom, other, keep_downstream_position = keep_downstream_position)
    end

    # During the last step, this atom might have been severed in an
    # inter-residue connection while being a child, therefore, it's parent was
    # assigned to the root. We should remove its own parent, therefore, in case
    # it happened (so it does not appear on root.children).
    _root = root(pose.graph)
    hasparent(atom) && isparent(_root, atom) && begin
        popparent!(atom)
        # hasparent(atom.container) && popparent!(atom.container)
    end

    # Remove from graph
    deleteat!(atom.container.items, findfirst(atom, atom.container.items))
    atom.container.size -= 1

    # Remove from state
    popped_state = splice!(pose.state, atom.index)

    # Reindex and set ascendents
    reindex(pose.graph) # Since we removed an atom, needs to be reindexed
    reindex(pose.state) # Since we removed an atom, needs to be reindexed
    ProtoSyn.request_i2c!(pose.state)

    # Update container 'itemsbyname'
    pop!(atom.container.itemsbyname, atom.name)

    # Set common ID
    popped_atom.id = popped_state.id = genid()

    return Pose(popped_atom, popped_state)
end

"""
    pop_residue!(pose::Pose{Topology}, residue::Residue; [keep_downstream_position::Bool = false])

Pop and return the desired [`Residue`](@ref) `residue` from the given
[`Pose`](@ref) `pose`. This is peformed by popping each [`Atom`](@ref) of the
[`Residue`](@ref) `residue` individually. If `keep_downstream_position` is set
to `true` (is, by default), the downstream [`Residue`](@ref) position is
maintained (by calling [`request_c2i!`](@ref) and [`sync!`](@ref) methods). In
either case, this method requests internal to cartesian coordinates conversion
at the end (using the [`request_i2c!`](@ref) method).

!!! ukw "Note:"
    The resulting [`Pose`](@ref) is re-indexed, therefore [`Residue`](@ref)
    `N + 1` becomes [`Residue`](@ref) `N`. When removing multiple
    [`Residue`](@ref) instances, consider performing a reversed loop.

# See also
[`pop_atom!`](@ref)

# Examples
```
julia> r = ProtoSyn.pop_residue!(pose, pose.graph[1][5])
Pose{Residue}(Residue{/ALA:51397}, State{Float64}:
 Size: 10
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function pop_residue!(pose::Pose{Topology}, residue::Residue; keep_downstream_position::Bool = true)

    if residue.container.container !== pose.graph
        error("Given Residue does not belong to the provided topology.")
    end

    # Instantiate the Residue and State to return
    popped_residue = Residue(residue.name.content, 1)
    popped_state   = State()

    # Copy residue to perserve parenthood and bonds information
    residue_copy = copy(residue)
    downstream = residue.children

    # Remove internal atoms. Notice the inverse loop. Also removes atom-level
    # parenthood and bonds
    for atom in reverse(residue.items)
        id = atom.id
        popped_atom = pop_atom!(pose, atom, keep_downstream_position = keep_downstream_position)
        insert!(popped_residue, 1, [popped_atom.graph]) # Save popped graph
        popped_residue[1].id = id
        append!(popped_state, popped_atom.state)        # Save popped state
    end

    # Recover intra-residue parenthoods and bonds in the popped_residue
    for (old_atom, new_atom) in zip(eachatom(residue_copy), eachatom(popped_residue))
        for child in old_atom.children
            child_id = child.id
            new_child = [atom for atom in eachatom(popped_residue) if atom.id == child_id][1]
            setparent!(new_child, new_atom)
        end

        for bond_atom in old_atom.bonds
            bond_id = bond_atom.id
            new_bond = [atom for atom in eachatom(popped_residue) if atom.id == bond_id][1]
            bond(new_atom, new_bond)
        end
    end

    # Remove from container.items
    deleteat!(residue.container.items, findfirst(residue, residue.container.items))
    residue.container.size -= 1

    # Remove parenthood to parent residue
    ProtoSyn.popparent!(residue)

    # Set downstream residue(s) parent to root
    _root = ProtoSyn.root(pose.graph).container
    for residue in downstream
        ProtoSyn.popparent!(residue)
        ProtoSyn.setparent!(residue, _root)
    end

    # Reindex
    reindex(pose.graph) # Since we removed a residue, needs to be reindexed.
    # Note that there's no need to reindex the state, it was reindexed in
    # pop_atom! and no state changes have occured since. The same occurs for
    # request_i2c! - this was performed in pop_atom! and no sync! call has
    # occured since, pose.state.i2c should be set to `true`.

    # Return the popped residue as a Fragment (no connection to root)
    popped_residue.id = popped_state.id = genid()
    return Pose(popped_residue, popped_state)
end


"""
# TODO
"""
function replace_by_fragment!(pose::Pose, atom::Atom, fragment::Fragment)
    
    # 0. Save current atom information
    fragment = copy(fragment)
    parent = atom.parent
    parent_index = parent.index
    parent_index_in_res = findfirst(x -> x === parent, parent.container.items)
    atomstate = copy(pose.state[atom])

    dihedrals = Vector{eltype(pose.state)}()
    _fragment = Pose(fragment)
    for a in ProtoSyn.root(_fragment.graph).children[1].children[1].children
        push!(dihedrals, ProtoSyn.getdihedral(_fragment.state, a))
    end
    
    # 1. Remove selected atom & any children atoms
    to_remove = ProtoSyn.travel_graph(atom)

    for a in reverse(to_remove) # Note the reverse loop
        ProtoSyn.pop_atom!(pose, a)
    end

    # 2. Add fragment to residue
    frag_origin = ProtoSyn.origin(fragment.graph)
    first_atom = frag_origin.children[1]
    insert!(parent.container, parent_index_in_res + 1, first_atom)
    ProtoSyn.popparent!(first_atom)
    ProtoSyn.unbond!(fragment, first_atom, frag_origin)
    ProtoSyn.bond(first_atom, parent)
    ProtoSyn.setparent!(first_atom, parent)
    atomstate.b = fragment.state[first_atom].b
    insert!(pose.state, parent_index + 1, State([atomstate]))

    for (i, frag_atom) in enumerate(ProtoSyn.travel_graph(first_atom)[2:end])
        insert!(parent.container, parent_index_in_res + i + 1, frag_atom)
        insert!(pose.state, parent_index + i + 1, State([fragment.state[frag_atom]]))
    end

    reindex(pose.graph, set_ascendents = true)
    reindex(pose.state)

    for (i, a) in enumerate(first_atom.children)
        ProtoSyn.setdihedral!(pose.state, a, dihedrals[i])
    end

    ProtoSyn.request_i2c!(pose.state)

    return pose
end


"""
# TODO
"""
function identify_atom_by_bonding_pattern(container::AbstractContainer, pattern::Vector{String})

    function seek_pattern(atom::Atom, inner_pattern::Vector{String}, previous::Opt{Atom})
        inner_candidates = [a for a in atom.bonds if (a.symbol === inner_pattern[1]) && (a !== previous)]
        length(inner_candidates) === 0 && return false
        deleteat!(inner_pattern, 1)
        length(inner_pattern) === 0 && return true
        return [seek_pattern(c, copy(inner_pattern), atom) for c in inner_candidates]
    end

    candidates = [a for a in eachatom(container) if a.symbol === pattern[1]]
    length(candidates) === 1 && return candidates[1]
    deleteat!(pattern, 1)
    length(pattern) === 0 && return candidates

    c = [seek_pattern(c, copy(pattern), nothing) for c in candidates]

    # Unravel results
    _c = Vector{Bool}([])
    for group in c
        for _ in 1:length(pattern)
            group = reduce(vcat, group)
        end
        push!(_c, any(group))
    end

    hits = candidates[_c]
    if length(hits) === 1
        return hits[1]
    else
        return hits
    end
end


"""
    sequence(container::ProtoSyn.AbstractContainer)::String
    sequence(pose::Pose)::String

Return the sequence of aminoacids (in 1 letter mode) of the given container/pose
as a string.

# Examples
```
julia> ProtoSyn.Peptides.sequence(pose)
"SESEAEFKQRLAAIKTRLQAL"
```
"""
function sequence(container::ProtoSyn.AbstractContainer)::String

    sequence = ""
    for residue in eachresidue(container)
        try
            sequence *= ProtoSyn.three_2_one[residue.name.content]
        catch KeyError
            sequence *= '?'
        end
    end

    return sequence
end

sequence(pose::Pose) = sequence(pose.graph)