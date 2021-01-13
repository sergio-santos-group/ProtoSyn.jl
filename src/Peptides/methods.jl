"""
    setss!(pose::Pose, (ϕ, ψ, ω)::NTuple{3,Number})

Set the `ϕ`, `ψ` and `ω` backbone angles of all residues in the given `pose`.
This function is usefull for setting the secondary structure of a pose. This
function acts on the internal coordinates and does not update cartesian
coordinates, although a request for conversion is made. It is up to the calling
function/user to explicitly synchornize coordinates via [`sync!`](@ref). In
certain cases it might be useful to not set the origin.children secondary
structure. An example is when a Segment is severed via [`unbond`](@ref), in
which case, updating the origin children will move one of the parts in an big
arm movement.
"""
function setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3, Number}, residues::Vector{Residue})
    state = container.state
    T = eltype(state)
    for r in residues
        if r.name == "PRO"
            # Proline is restricted to TRANS conformation
            # This conformation is most abundant (~95%) in globular proteins.
            # https://pubs.acs.org/doi/10.1021/jacs.0c02263
            # PHI: -75° | PSI: 145° | OMEGA: 180°
            ProtoSyn.setdihedral!(container.state, Dihedral.phi(r), T(-1.308997))
            # Last residues of chain might not have a psi angle.
            Dihedral.psi(r) !== nothing && ProtoSyn.setdihedral!(container.state, Dihedral.psi(r),  T(2.5307274))
            ProtoSyn.setdihedral!(container.state, Dihedral.omega(r), T(3.1415927))
            continue
        end
        ProtoSyn.setdihedral!(container.state, Dihedral.phi(r), ϕ)
        # Last residues of chain might not have a psi angle.
        Dihedral.psi(r) !== nothing && ProtoSyn.setdihedral!(container.state, Dihedral.psi(r),  ψ)
        ProtoSyn.setdihedral!(container.state, Dihedral.omega(r), ω)
    end

    return container
end

setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3, Number}) = begin
    residues::Vector{Residue} = collect(eachresidue(container.graph))
    Peptides.setss!(container, (ϕ, ψ, ω), residues)
end

setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3,Number}, sele::ProtoSyn.AbstractSelection) = begin
    residues = ProtoSyn.PromoteSelection(sele, Residue, any)(container, gather = true)
    Peptides.setss!(container, (ϕ, ψ, ω), residues)
end

setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3,Number}, residue::Residue) = begin
    Peptides.setss!(container, (ϕ, ψ, ω), [residue])
end


"""
    append_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ss::NTuple{3,Number} = SecondaryStructure[:linear], op = "α")

Based on the provided `grammar`, add the residue sequence from `derivation` to
the given `pose`, appending it AFTER the given `residue`. This residue and the
new fragment will be connected using operation `op` ("α" by default). Set the
secondary structure of the added residue to match the given `ss` (linear, by
default). Return the altered `pose`.

# Examples
```jldoctest
julia> append_residues!(pose, pose.graph[1][1], reslib, seq"A")
```
"""
function append_residues!(pose::Pose{Topology}, residue::Residue,
    grammar::LGrammar, derivation;
    ss::NTuple{3,Number} = SecondaryStructure[:linear], op = "α")

    Builder.append_residues!(pose, residue, grammar, derivation; op = op)
    residues = residue.container.items[(residue.index + 1):(residue.index + length(derivation))]
    setss!(pose, ss, residues)
    return pose
end


"""
    insert_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ss::NTuple{3, Number} = SecondaryStructure[:linear], op = "α")

Based on the provided `grammar`, add the residue sequence from `derivation` to
the given `pose`, inserting it ON THE POSITION of the given `residue` (the
`residue` gets shifted downstream). The first downstream residue and the new
fragment will be connected using operation `op` ("α" by default). Upstream
residues are also connected using this operation if they are not origin. Set the
secondary structure of the added residue to match the given `ss` (linear, by
default). Return the altered `pose`.

# Examples
```jldoctest
julia> insert_residues!(pose, pose.graph[1][2], reslib, seq"A")
```
"""
function insert_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation;
    ss::NTuple{3, Number} = SecondaryStructure[:linear], op = "α")

    connected_to_origin = residue.parent == ProtoSyn.origin(pose.graph).container
    if connected_to_origin
        state = pose.state
        N = state[residue["N"]] # This is the N atom state
        i = residue.index       # This is the residue index
        (b, θ, ϕ) = (N.b, N.θ, N.ϕ)
    else
        ProtoSyn.unbond(pose, residue.container[residue.index - 1]["C"], residue["N"])
    end
    
    Builder.insert_residues!(pose, residue, grammar, derivation; op = op,
        connect_upstream = !connected_to_origin)

    if connected_to_origin
        N = state[residue.container[i]["N"]]
        (N.b, N.θ, N.ϕ) = (b, θ, ϕ)
        ProtoSyn.request_i2c(state; all = true)
    end

    # residues = residue.container.items[(residue.index - length(derivation)):(residue.index)]
    setss!(pose, ss, residue.parent)
end


"""
    force_mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation, op = "α")

Mutate the given `pose` at `residue`, changing it's aminoacid to be
`derivation`, as given by the template at `grammar`. This function changes the
whole residue (backbone included). By default, the user should use `mutate!`
instead of this function, except for uncommon aminoacids.

# Examples
```jldoctest
julia> Peptides.force_mutate!(pose, pose.graph[1][3], res_lib, seq"K")
```
"""
function force_mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation, op = "α")

    # 1) Insert new residue in position R, shifting all downstream residues + 1
    # The inserted residue backbone dihedrals should match the secondary
    # structure pre-existent in the pose
    phi   = ProtoSyn.getdihedral(pose.state, Peptides.Dihedral.phi(residue))
    psi   = ProtoSyn.getdihedral(pose.state, Peptides.Dihedral.psi(residue))
    omega = ProtoSyn.getdihedral(pose.state, Peptides.Dihedral.omega(residue))
    old_r_pos = findfirst(residue, residue.container.items)
    Peptides.insert_residues!(pose, residue, grammar, derivation,
        ss = (phi, psi, omega), op = op)
    new_residue = residue.container.items[old_r_pos]

    # 2) Remove old residue (now in position R + 1)
    old_r_pos += 1
    Peptides.pop_residue!(pose, residue.container.items[old_r_pos])

    # 3) Remove downstream residue parents (now in position R + 1, since we
    # deleted the R + 1 old residue)
    downstream_r_pos = old_r_pos
    new_residue.container.items[downstream_r_pos].parent = nothing
    new_residue.container.items[downstream_r_pos]["N"].parent = nothing

    # 4) Use operator 'op' from 'grammar' to set parenthoods, bonds and correct
    # position, based on the peptidic bond
    grammar.operators[op](new_residue, pose, residue_index = downstream_r_pos)

    # 5) Re-set ascedents
    ProtoSyn.reindex(pose.graph)
    return pose
end


"""
    mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation)

Mutate the given `pose` at `residue`, changing it's aminoacid to be
`derivation`, as given by the template at `grammar`. This function changes the
sidechain only (± 7x faster than `force_mutate!`). When mutating to Proline,
falls back to `force_mutate!`. *Note:* sidechains are selected based on the atom
name (backbone atoms must be named N, H, CA, C and O, exclusively. Non backbone
atoms should have other names, such as H1, H2, etc.)

# Examples
```jldoctest
julia> Peptides.mutate!(pose, pose.graph[1][3], res_lib, seq"K")
```
"""
function mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation)

    @assert length(derivation) == 1 "Derivation must have length = 1."
    
    sidechain = (!an"^CA$|^N$|^C$|^H$|^O$"r)(residue, gather = true)

    same_aminoacid = string(Peptides.three_2_one[residue.name]) == derivation[1]
    if same_aminoacid && length(sidechain) > 0
        println("No mutation required, residue already has sidechain of the requested type.")
        return pose
    end

    if derivation[1] == "P"
        return Peptides.force_mutate!(pose, residue, grammar, derivation)
    end

    frag = Builder.fragment(grammar, derivation)
    
    # Remove old sidechain
    for atom in sidechain
        Builder.pop_atom!(pose, atom)
    end
    
    # Insert new sidechain
    frag_sidechain = (!an"^CA$|^N$|^C$|^H$|^O$"r)(frag, gather = true)
    poseCA = residue["CA"]
    
    # poseCA_index is the LOCAL index (inside residue.items). poseCA.index is
    # the index in the whole pose
    poseCA_index = findfirst(x -> x === poseCA, residue.items)
    
    objective_changes = []
    # The ϕ is the Phi of the fragment, reduced to be in [0, 360] degrees range.
    # All CA child atom's dihedrals will be placed in relationship (relative) to
    # the corresponding residue Phi dihedral.
    # Here we are measuring the default difference between those two angles in
    # the template fragment. This value could, in a later version of ProtoSyn,
    # be parametrized somewhere.
    _ϕ = ProtoSyn.getdihedral(frag.state, Peptides.Dihedral.phi(frag.graph[1]))
    ϕ  = ProtoSyn.unit_circle(_ϕ)
    for (index, atom) in enumerate(frag_sidechain)
        parent_is_CA = false
        if atom.parent.name == "CA"
            parent_is_CA = true
            _atom_dihedral = ProtoSyn.getdihedral(frag.state, atom)
            atom_dihedral  = ProtoSyn.unit_circle(_atom_dihedral)
            objective_change = atom_dihedral - ϕ
            push!(objective_changes, objective_change)
            ProtoSyn.unbond(frag, atom, atom.parent)
        end

        # Insert into the graph
        # Note: insert! already sets the residue.itemsbyname
        insert!(residue, poseCA_index + index, atom)

        # Since now atom is already in the graph, we can bond and add parents
        if parent_is_CA
            ProtoSyn.bond(atom, poseCA)
            setparent!(atom, poseCA)
        end
    end

    _start = frag_sidechain[1].index
    _end   = frag_sidechain[end].index
    insert!(pose.state, poseCA.index + 1, splice!(frag.state, _start:_end))
    
    reindex(pose.graph)

    # Fix CA children positions
    pose_sidechain = (!an"^CA$|^N$|^C$|^H$|^O$"r)(residue, gather = true)
    Δϕ             = pose.state[residue["CA"]].Δϕ
    index          = 1
    _ϕ = ProtoSyn.getdihedral(pose.state, Peptides.Dihedral.phi(residue))
    ϕ  = ProtoSyn.unit_circle(_ϕ)
    for child in residue["CA"].children
        if child in pose_sidechain
            objective = ϕ + objective_changes[index]
            pose.state[child].ϕ = ProtoSyn.unit_circle(objective - Δϕ)
            index += 1
        end
    end
    
    residue.name = Peptides.one_2_three[derivation[1][1]]
    ProtoSyn.request_i2c(pose.state)

    return pose
end


"""
    pop_residue!(pose::Pose{Topology}, residue::Residue)

Delete `residue` from the given `pose`. Removes parenthood (at residue and atom
level), removes bonds, reindexes ascedents and sets position and parenthood to
origin of downstream residues.

# Examples
```jldoctest
julia> Peptides.pop_residue!(pose, pose.graph[1][3])
```
"""
function pop_residue!(pose::Pose{Topology}, residue::Residue)

    # 1) Since we are calling Peptides, we know that we should unbond this
    # residue and the next (children)
    for child in residue.children
        Peptides.unbond(pose, residue, child) # Order is important
    end

    # 2) Now we can safelly pop the residue, while maintaining the positions of
    # downstream residues
    Builder.pop_residue!(pose, residue)
end


"""
    remove_sidechains!(pose::Pose{Topology}, selection::Opt{AbstractSelection} = nothing)

Removes the sidechain atoms of the given `pose`. If a `selection` is given, only
the atom of that selection are considered for possible removal.

# Examples
```jldoctest
julia> Peptides.remove_sidechains!(pose)

julia> Peptides.remove_sidechains!(pose, rn"ALA")
```
"""
function remove_sidechains!(pose::Pose{Topology}, selection::Opt{AbstractSelection} = nothing)
    _selection = !(an"^CA$|^N$|^C$|^H$|^O$"r | rn"PRO")
    if selection !== nothing
        _selection = _selection & selection
    end
    sidechain = _selection(pose, gather = true)
    for atom in reverse(sidechain) # Note the reverse loop
        Builder.pop_atom!(pose, atom)
    end

    return pose
end


"""
    add_sidechains!(pose::Pose{Topology}, grammar::LGrammar, selection::Opt{AbstractSelection} = nothing)

Add the sidechain atoms to the given `pose`, based on the templates of
`grammar`. If a `selection` is given, only the residues of that selection
(promoted to `Residue` instances, using the default aggregator function) are
considered for sidechain addition. The addition is performed using `mutate!`
function.

# Examples
```jldoctest
julia> Peptides.add_sidechains!(pose)

julia> Peptides.add_sidechains!(pose, rn"ALA")
```
"""
function add_sidechains!(pose::Pose{Topology}, grammar::LGrammar, selection::Opt{AbstractSelection} = nothing)
    if selection !== nothing
        residues = selection(pose, gather = true)
    else
        residues = collect(eachresidue(pose.graph))
    end
    for residue in residues
        derivation = [string(Peptides.three_2_one[residue.name])]
        Peptides.mutate!(pose, residue, grammar, derivation)
    end

    return pose
end


"""
    unbond(pose::Pose, residue_1::Residue, residue_2::Residue)

Unbond the two provided residues. In order to do this, perform the following
steps, as in ProtoSyn.unbond:
 - Unset parents/children
 - Unbond neighbours
 - Remove from graph
 - Remove from state
 - Set new ascendents
 - Update the container itemsbyname
Furthermore, since the peptidic nature is known:
 - Set correct positioning of downstream residues

# Examples
```jldoctest
julia> unbond(pose, pose.graph[1][2])
```
"""
function unbond(pose::Pose, residue_1::Residue, residue_2::Residue)

    isparent(residue_1, residue_2) && return _unbond(pose, residue_1, residue_2)
    isparent(residue_2, residue_1) && return _unbond(pose, residue_2, residue_1)
end

function _unbond(pose::Pose, residue_1::Residue, residue_2::Residue)

    ProtoSyn.unbond(pose, residue_1["C"], residue_2["N"])
    
    # Set correct positioning
    state = pose.state
    _origin = ProtoSyn.origin(pose.graph)
    sync!(pose)

    at_N = state[residue_2["N"]]
    at_N.b = ProtoSyn.distance(at_N, state[_origin])
    at_N.θ = ProtoSyn.angle(at_N, state[_origin], state[_origin.parent])
    v = ProtoSyn.dihedral(at_N, state[_origin], state[_origin.parent], state[_origin.parent.parent])
    at_N.ϕ += v - ProtoSyn.getdihedral(state, residue_2["N"])

    for child in residue_2["N"].children
        at_1 = state[child]
        at_1.θ = ProtoSyn.angle(at_1, at_N, state[_origin])
        v = ProtoSyn.dihedral(at_1, at_N, state[_origin], state[_origin.parent])
        at_1.ϕ += v - ProtoSyn.getdihedral(state, child)

        for grandchild in child.children
            at_2 = state[grandchild]
            v = ProtoSyn.dihedral(at_2, at_1, at_N, state[_origin])
            at_2.ϕ += v - ProtoSyn.getdihedral(state, grandchild)
        end
    end

    ProtoSyn.request_i2c(state, all = true)
end