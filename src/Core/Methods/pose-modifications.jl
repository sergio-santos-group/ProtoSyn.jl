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
    add_hydrogens!(pose::Pose, res_lib::LGrammar, selection::Opt{AbstractSelection} = nothing)

Predict and add hydrogens to the provided [`Pose`](@ref) `pose`, according to
the templates in the given [`LGrammar`](@ref) `res_lib`. Note that both residues
and atoms are retrieved from the `res_lib` based on their name. If provided,
an `AbstractSelection` `selection` limits the selected atoms to receive
hydrogens (defaults to an [`Atom`](@ref) level `AbstractSelection`, expects
[`Atom`](@ref) instances in the given [`Pose`](@ref) `pose` to be correctly
indexed). This function applies internal coordinates based on the template from
`res_lib` where the dihedral angle is rotated in order to match the current
dihedral angles in the [`Pose`](@ref) `pose`. Make sure the internal coordinates
are synched (using the [`sync!`](@ref) method).

# Examples
```jldoctest
julia> ProtoSyn.add_hydrogens!(pose, Peptides.grammar, rid"1:10")
Pose{Topology}(Topology{/1bkr:47568}, State{Float64}:
 Size: 975
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function add_hydrogens!(pose::Pose, res_lib::LGrammar, selection::Opt{AbstractSelection} = nothing)

    if selection !== nothing
        atoms = ProtoSyn.promote(selection, Atom)(pose, gather = true)
    else
        atoms = TrueSelection{Atom}()(pose, gather = true)
    end

    # Gather all atom names involved in inter-residue bonds (from the LGrammar
    # operators list)
    r1s = [res_lib.operators[op].r1 for op in keys(res_lib.operators)]
    r2s = [res_lib.operators[op].r2 for op in keys(res_lib.operators)]

    i = 1 # For hydrogen naming. It's recommender to re-name atoms after
    for residue in eachresidue(pose.graph)

        # Create the residue template from residue name (with full hydrogens)
        code = string(ProtoSyn.three_2_one[residue.name])
        template = ProtoSyn.getvar(res_lib, code)

        for atom in residue.items
            
            # Ignore hydrogen already placed in the original pose
            atom.symbol === "H" && continue

            # Ignore non-selected atoms
            !(atom in atoms) && continue
            
            # Get the corresponding template atom (by name) 
            temp                = template.graph[1][atom.name]

            # Predict the current coordination (number of bonds) for this atom.
            # 1. Check inter-residue operators in the provided LGrammar to see
            # if this atom (in the template) can be involved in any
            # inter-residue bond. This doesn't mean that the corresponding atom
            # in the original pose is indeed involved in an inter-residue bond,
            # but the current version of ProtoSyn assumes that it is. In a
            # future iteration, checking the atom-pair distance to atoms in
            # other residues (with the right name, as defined in the LGrammar
            # operation) might help to better determine inter-residue bonds.
            inter_residue_bonds = 0
            if atom.name in r1s
                inter_residue_bonds += 1
            end
            if atom.name in r2s
                inter_residue_bonds += 1
            end
            # 2. Add the intra-residue bonds. The number of hydrogens missing
            # (N) is defined as the difference between the number of bonds in
            # the template and the original pose (considering that all atoms
            # that can form inter-residue bonds do).
            temp_bonds = length(temp.bonds) + inter_residue_bonds
            N          = temp_bonds - length(atom.bonds)

            # If no hydrogens are missing, ignore next steps.
            N < 1 && continue

            # The missing hydrogens will be copied from the template (via the
            # internal coordinates). However, the internal coordinates in the
            # template are only in accordance to the template dihedral angles.
            # Any different dihedral angles (phi, psi, omega and chis) wil break
            # the structure and induce super-positions. The final internal
            # coordinates to the introduced hydrogens need to reflect the
            # current dihedral angles the parent atom is involved.
            # 1. Measure current dihedral angle on the parent atom. This depends
            # on the chosen child atom to measure the dihedral. Since the final
            # internal coordinates of the newly introduced hdyrogens will be
            # "relative" to this measured dihedral, it is indiferent which child
            # atom is picked for measurement. Note that this function expects
            # the original pose to be synched (will measure dihedrals from
            # internal coordinates). If the currently focused atom doesn't have
            # children, then the dihedral angle doesn't have any interference
            # with the default tempalte internal coordinates (the relative
            # difference to apply is 0.0).
            if length(atom.children) === 0
                Δχ = 0.0
            else
                pose_child = atom.children[1]
                χ_pose     = ProtoSyn.getdihedral(pose.state, pose_child)
                temp_child = temp.children[findfirst((a) -> a.name === pose_child.name, temp.children)]
                χ_temp     = ProtoSyn.getdihedral(template.state, temp_child)
                Δχ         = χ_pose - χ_temp
            end

            for t in temp.bonds
                t.symbol !== "H" && continue

                temp_state = template.state[t]

                h = Atom("H$i", -1, -1, "H")
                h_state = AtomState()
                h_state.b = temp_state.b
                h_state.θ = temp_state.θ
                h_state.ϕ = temp_state.ϕ - pose.state[atom.parent].Δϕ
                h_state.ϕ = h_state.ϕ + Δχ
                ProtoSyn.insert_atom_as_children!(pose, atom, h, h_state)
                
                i += 1
            end
        end
    end

    reindex(pose.graph; set_ascendents = true)
    reindex(pose.state)
    ProtoSyn.request_i2c!(pose.state; all = true)
    sync!(pose)
end


"""
    replace_by_fragment!(pose::Pose, atom::Atom, fragment::Fragment)

Replace the selected [`Atom`](@ref) `atom` instance (and any downstream children
atom, uses [`travel_graph`](@ref) stating on [`Atom`](@ref) `atom` to define the
downstream region) with the given [`Fragment`](@ref) `fragment`, in the context
of the provided [`Pose`](@ref) `pose` (updates the [`State`](@ref) and
[Graph](@ref), etc). The first [`Atom`](@ref) in the [`Fragment`](@ref)
`fragment` (also known as root or R [`Atom`](@ref)) is placed in the same
position as the chosen [`Atom`](@ref) `atom` for replacement, and then removed.
This serves only to orient the remaining [`Fragment`](@ref) `fragment`. Requests
internal to cartesian coordinates.

# Examples
```
julia> frag = ProtoSyn.getvar(ProtoSyn.modification_grammar, "PO4")
Fragment(Segment{/po4:18413}, State{Float64}:
 Size: 5
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)

julia> ProtoSyn.replace_by_fragment!(pose, pose.graph[1, 3, "HG"], frag)
Pose{Topology}(Topology{/2a3d:3900}, State{Float64}:
 Size: 1143
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
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
    first_atom  = frag_origin.children[1]
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
    if atom.name in keys(atom.container.itemsbyname)
        pop!(atom.container.itemsbyname, atom.name)
    else
        @warn "While popping atom $(atom.name), attempted to remove entry from the container.itemsbyname, but no entry was found."
    end

    # Set common ID
    popped_atom.id = popped_state.id = genid()

    return Pose(popped_atom, popped_state)
end


"""
    pop_atoms!(pose::Pose{Topology}, selection::Opt{AbstractSelection} = nothing; keep_downstream_position::Bool = false)

Remove all selected atoms (by the given `AbstractSelection` `selection`) from
the provided [`Pose`](@ref) `pose`. `keep_downstream_position` sets whether to
re-calculate the internal coordinates of downstream [`Atom`](@ref) instances
(children) from cartesian coordinates (set to `false`, by default).

# See also
[`pop_atom!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.pop_atoms!(pose, as"H")
Pose{Topology}(Topology{/1bkr:24417}, State{Float64}:
 Size: 887
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function pop_atoms!(pose::Pose{Topology}, selection::Opt{AbstractSelection} = nothing; keep_downstream_position::Bool = false)
    
    if selection !== nothing
        atoms = ProtoSyn.promote(selection, Atom)(pose, gather = true)
    else
        atoms = TrueSelection{Atom}(pose, gather = true)
    end

    for atom in reverse(atoms)
        ProtoSyn.pop_atom!(pose, atom; keep_downstream_position = keep_downstream_position)
    end

    return pose
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