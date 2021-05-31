using ProtoSyn.Units

# * ----------------------------------------------------------------------------
# * Append and insert fragments

"""
    append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, frag::Pose{Segment}; [ss::Opt{SecondaryStructureTemplate} = nothing], [op = "α"])

Add the [`Fragment`](@ref) `frag` to the given [`Pose`](@ref) `pose`, appending
it after the given [`Residue`](@ref) `residue`. This residue and the new
[`Fragment`](@ref) `frag` will be connected using operation `op` ("α" by
default) of the given [`LGrammar`](@ref) `grammar`. If given, a
[`SecondaryStructureTemplate`](@ref) `ss` can be applied to the new appendage
(using the [`setss!`](@ref) method). In either case, the `C=O` bond position is
re-calculated and set (in the first residue of the appendage). Request internal
to cartesian coordinate conversion and return the altered [`Pose`](@ref) `pose`.

!!! ukw "Note:"
    This function is an extension of [`ProtoSyn.append_fragment!`](@ref).

# See also
[`insert_fragment!`](@ref ProtoSyn.Peptides.insert_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Pose{Segment}; ::Opt{SecondaryStructureTemplate}, ::Any))

# Examples
```jldoctest
julia> ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][end], res_lib, frag)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 373
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, frag::Pose{Segment}; ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")
    
    ProtoSyn.append_fragment!(pose, residue, grammar, frag, op = op)
    sync!(pose) # Apply changes from `append_fragment!`

    # By default, the grammar operator alpha sets the psi dihedral (between the
    # anchor aminoacid and the first aminoacid of the fragment) to 180°. This
    # usually interfered with the current position of the oxygen atom (in the
    # C=O bond. For that reason, either the oxygen atom or the nitrogen atom
    # need to be rotated. In the current version of this method, the oxygen atom
    # position is rotated to accomodate the fragment appendage, by setting its
    # ϕ value to be 180° rotated from the current psi dihedral.
    current_psi = ProtoSyn.getdihedral(pose.state, residue.children[1]["N"])
    target_dihedral = current_psi - 180°
    pose.state[residue["O"]].ϕ = pose.state[residue["C"]].Δϕ - target_dihedral
    ProtoSyn.request_i2c!(pose.state, all = true)
    
    if ss !== nothing
        residues = residue.container.items[(residue.index + 1):(residue.index + length(frag.graph))]
        ProtoSyn.Peptides.setss!(pose.state, ss, residues)
    end

    return pose
end


"""
    append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; [ss::Opt{SecondaryStructureTemplate} = nothing], [op = "α"])

Add the a [`Fragment`](@ref) built from the given `derivation` to the provided
[`Pose`](@ref) `pose`, appending it after the given [`Residue`](@ref) `residue`.
This residue and the new [`Fragment`](@ref) `frag` will be connected using
operation `op` ("α" by default) of the given [`LGrammar`](@ref) `grammar`. If
given, a [`SecondaryStructureTemplate`](@ref) `ss` can be applied to the new
appendage (using the [`setss!`](@ref) method). In either case, the `C=O` bond
position is re-calculated and set (in the anchor for the first residue of the
appendage). Request internal to cartesian coordinate conversion and return the
altered [`Pose`](@ref) `pose`.

!!! ukw "Note:"
    This function is an extension of [`ProtoSyn.append_fragment!`](@ref).

# See also
[`insert_fragment!`](@ref ProtoSyn.Peptides.insert_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Pose{Segment}; ::Opt{SecondaryStructureTemplate}, ::Any))

# Examples
```jldoctest
julia> ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][end], res_lib, seq"AAA")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 373
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")
    frag = ProtoSyn.fragment(grammar, derivation)
    return ProtoSyn.Peptides.append_fragment!(pose, residue, grammar, frag, ss = ss,
        op = op)
end


"""
    insert_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, frag::Pose{Segment}; ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")

Insert the [`Fragment`](@ref) `frag` in the given `pose`, on the position of the
provided [`Residue`](@ref) instance `residue` (the `residue` gets shifted
downstream). This first downstream [`Residue`](@ref) and the new
[`Fragment`](@ref) will be connected using operation `op` ("α" by default) from
[`LGrammar`] `grammar`. Also connects to the upstream [`Residue`](@ref)
instance, using the same operation. If given, a
[`SecondaryStructureTemplate`](@ref) `ss` can be applied to the new appendage
(using the [`setss!`](@ref) method). If the appendage is not being inserted at
the [`root`](@ref), the `C=O` bond position is re-calculated and set (in the
anchor for the first residue of the appendage). If the appendage is being
inserted at the [`root`](@ref), perform a soft uncap of the terminal hydrogen
atoms (removes "H2" and "H3", leaves "H1", renames it to "H") and recalculate
the N-H bond position (at the first downstream [`Residue`](@ref)). Request
internal to cartesian coordinate conversion and return the altered
[`Pose`](@ref) `pose`.

!!! ukw "Note:"
    This function is an extension of [`ProtoSyn.insert_fragment!`](@ref).

# See also
[`append_fragment!`](@ref ProtoSyn.append_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Pose{Segment}; ::Opt{SecondaryStructureTemplate}, ::Any))

# Examples
```jldoctest
julia> ProtoSyn.Peptides.insert_fragment!(pose, pose.graph[1][1], res_lib, frag)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 373
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function insert_fragment!(pose::Pose{Topology}, residue::Residue,
    grammar::LGrammar, frag::Pose{Segment};
    ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")

    # Setup
    frag_size = length(eachresidue(frag.graph))
    anchor    = residue.parent
    r_id      = residue.index

    # Case being connected to root, perform soft uncap the N-terminal
    anchor == ProtoSyn.root(pose.graph).container && begin
        residue["H2"] !== nothing && ProtoSyn.pop_atom!(pose, residue["H2"])
        residue["H3"] !== nothing && ProtoSyn.pop_atom!(pose, residue["H3"])

        # Rename H1 to H
        if residue["H1"] !== nothing
            residue.itemsbyname["H"] = residue["H1"]
            pop!(residue.itemsbyname, "H1")
            residue["H"].name = "H"
        end
    end

    # Insert the fragments
    ProtoSyn.insert_fragment!(pose, residue, grammar, frag, op = op)

    if !(anchor == ProtoSyn.root(pose.graph).container)
        # If not connected to root, adjust the position of the upstream C=O bond
        # current_psi = ProtoSyn.getdihedral(pose.state, anchor.children[1]["N"])
        # target_dihedral = current_psi - 180°
        # pose.state[anchor["O"]].ϕ = pose.state[anchor["C"]].Δϕ - target_dihedral
        # ProtoSyn.request_i2c!(pose.state, all = true)
    else
        # If connected to root, adjust the position of the downstream N=H bond
        current_omega = ProtoSyn.getdihedral(pose.state, residue["CA"])
        target_dihedral = current_omega - 180°
        pose.state[residue["H"]].ϕ = pose.state[residue["N"]].Δϕ - target_dihedral
        ProtoSyn.request_i2c!(pose.state, all = true)
    end

    # Apply secondary structure
    if ss !== nothing
        residues = residue.container.items[r_id:r_id + frag_size - 1]
        ProtoSyn.Peptides.setss!(pose.state, ss, residues)
    end

    return pose
end


"""
    insert_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")

Insert the [`Fragment`](@ref) built from the given `derivation` in the provided
`pose`, on the position of the provided [`Residue`](@ref) instance `residue`
(the `residue` gets shifted downstream). This first downstream [`Residue`](@ref)
and the new [`Fragment`](@ref) will be connected using operation `op` ("α" by
default) from [`LGrammar`] `grammar`. Also connects to the upstream
[`Residue`](@ref) instance, using the same operation. If given, a
[`SecondaryStructureTemplate`](@ref) `ss` can be applied to the new appendage
(using the [`setss!`](@ref) method). If the appendage is not being inserted at
the [`root`](@ref), the `C=O` bond position is re-calculated and set (in the
anchor for the first residue of the appendage). If the appendage is being
inserted at the [`root`](@ref), perform a soft uncap of the terminal hydrogen
atoms (removes "H2" and "H3", leaves "H1", renames it to "H") and recalculate
the N-H bond position (at the first downstream [`Residue`](@ref)). Request
internal to cartesian coordinate conversion and return the altered
[`Pose`](@ref) `pose`.

!!! ukw "Note:"
    This function is an extension of [`ProtoSyn.insert_fragment!`](@ref).

# See also
[`append_fragment!`](@ref ProtoSyn.append_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Pose{Segment}; ::Opt{SecondaryStructureTemplate}, ::Any))

# Examples
```jldoctest
julia> ProtoSyn.Peptides.insert_fragment!(pose, pose.graph[1][1], res_lib, seq"AAA")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 373
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function insert_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")
    frag = ProtoSyn.fragment(grammar, derivation)
    return ProtoSyn.Peptides.insert_fragment!(pose, residue, grammar, frag, ss = ss,
        op = op)
end

# * ----------------------------------------------------------------------------
# * Mutate aminoacids

"""
    mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation)

Mutate the given [`Pose`](@ref) `pose` at [`Residue`](@ref) `residue`, changing
it's aminoacid to be `derivation`, as given by the template at `grammar`. This
function changes the sidechain only (± 7x faster than [`force_mutate!`](@ref)).
When mutating to Proline, falls back to [`force_mutate!`](@ref). If
`ignore_existing_sidechain` is set to `true` (`false` by default), existing
sidechains are first removed and then re-added, regardless of being of the same
type (which normally are ignored, if no mutation is required).

!!! ukw "Note:"
    Sidechains are selected based on the [`Atom`](@ref)`.name` (backbone
    [`Atom`](@ref) instances must be named N, H, CA, C and O, exclusively. Non
    backbone [`Atom`](@ref) instances should have other names, such as H1, H2,
    etc.)

# Examples
```jldoctest
julia> ProtoSyn.Peptides.mutate!(pose, pose.graph[1][3], res_lib, seq"K")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 354
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ignore_existing_sidechain::Bool = false)

    @assert length(derivation) == 1 "Derivation must have length = 1."
    
    sidechain = (!an"^CA$|^N$|^C$|^H$|^O$"r)(residue, gather = true)

    same_aminoacid = string(ProtoSyn.Peptides.three_2_one[residue.name]) == derivation[1]
    if same_aminoacid && length(sidechain) > 0 && !ignore_existing_sidechain
        # println("No mutation required, residue already has sidechain of the requested type.")
        return pose
    end

    if derivation[1] == "P"
        return Peptides.force_mutate!(pose, residue, grammar, derivation)
    end

    frag = ProtoSyn.fragment(grammar, derivation)
    
    # Remove old sidechain
    for atom in reverse(sidechain)
        ProtoSyn.pop_atom!(pose, atom; keep_downstream_position = false)
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
            parent_is_CA     = true
            _atom_dihedral   = ProtoSyn.getdihedral(frag.state, atom)
            atom_dihedral    = ProtoSyn.unit_circle(_atom_dihedral)
            objective_change = atom_dihedral - ϕ
            push!(objective_changes, objective_change)
            ProtoSyn.unbond!(frag, atom, atom.parent)
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
    reindex(pose.state)

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
    
    residue.name = ProtoSyn.ResidueName(Peptides.one_2_three[derivation[1][1]])
    ProtoSyn.request_i2c!(pose.state)

    return pose
end


"""
    force_mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation, op = "α")

Mutate the given [`Pose`](@ref) `pose` at [`Residue`](@ref) `residue`, changing
it's aminoacid to be `derivation`, as given by the template at `grammar`. This
function changes the whole residue (backbone included). By default, the user
should use [`mutate!`](@ref) instead of this function, except for uncommon
aminoacids.

# Examples
```jldoctest
julia> ProtoSyn.Peptides.force_mutate!(pose, pose.graph[1][3], res_lib, seq"K")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 354
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function force_mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation, op = "α")

    # 1) Insert new residue in position R, shifting all downstream residues + 1
    # The inserted residue backbone dihedrals should match the secondary
    # structure pre-existent in the pose
    phi   = ProtoSyn.getdihedral(pose.state, Peptides.Dihedral.phi(residue))
    psi   = ProtoSyn.getdihedral(pose.state, Peptides.Dihedral.psi(residue))
    omega = ProtoSyn.getdihedral(pose.state, Peptides.Dihedral.omega(residue))
    ss    = ProtoSyn.Peptides.SecondaryStructureTemplate(phi, psi, omega)
    old_r_pos = findfirst(residue, residue.container.items)
    Peptides.insert_fragment!(pose, residue, grammar, derivation,
        ss = ss, op = op)
    new_residue = residue.container.items[old_r_pos]

    # 2) Remove old residue (now in position R + 1)
    old_r_pos += 1
    # N = residue.container.items[old_r_pos]["N"]
    ProtoSyn.pop_residue!(pose, residue.container.items[old_r_pos],
        keep_downstream_position = false)
    # i = findall(a -> a == N, ProtoSyn.root(pose.graph).children)
    # deleteat!(ProtoSyn.root(pose.graph).children, i)

    # 3) Remove downstream residue parents (now in position R + 1, since we
    # deleted the R + 1 old residue). Current parents (after pop_residue! are
    # the root).
    downstream_r_pos = old_r_pos
    ProtoSyn.popparent!(new_residue.container.items[downstream_r_pos])
    ProtoSyn.popparent!(new_residue.container.items[downstream_r_pos]["N"])

    # 4) Use operator 'op' from 'grammar' to set parenthoods, bonds and correct
    # position, based on the peptidic bond
    grammar.operators[op](new_residue, pose, residue_index = downstream_r_pos)

    # 5) Re-set ascedents
    ProtoSyn.reindex(pose.graph)
    ProtoSyn.request_i2c!(pose.state, all = true)

    return pose
end