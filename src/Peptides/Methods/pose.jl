using ProtoSyn.Units
using ProgressMeter

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
    elseif residue["H"] !== nothing
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

    mutate!(pose::Pose{Topology}, sele::AbstractSelection, grammar::LGrammar, derivation)

Mutate the given [`Pose`](@ref) `pose` at [`Residue`](@ref) `residue` selected
by the `AbstractSelection` `sele`. Note that this selection can only return a
single aminoacid at a time.

!!! ukw "Note:"
    Sidechains are selected based on the [`Atom`](@ref)`.name` (backbone
    [`Atom`](@ref) instances must be named N, H, CA, C and O, exclusively. Non
    backbone [`Atom`](@ref) instances should have other names, such as HB1, HB2,
    etc.) This function attempts to keep terminal caps intact (doesn't remove
    any atom named H1, H2, etc or OXT).

# See also
[`force_mutate!`](@ref)

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

    T = eltype(pose.state)
    @assert length(derivation) == 1 "Derivation must have length = 1."
    
    # ! The definition of sidechain could be improved in future versions !
    sidechain = (!an"^CA$|^N$|^C$|^H\d*$|^O$|^OXT$"r)(residue, gather = true)

    same_aminoacid = string(ProtoSyn.ProtoSyn.three_2_one[residue.name]) == derivation[1]
    if same_aminoacid && length(sidechain) > 0 && !ignore_existing_sidechain
        # println("No mutation required, residue already has sidechain of the requested type.")
        return pose
    end

    if derivation[1] == "P"
        return Peptides.force_mutate!(pose, residue, grammar, derivation)
    end

    frag = ProtoSyn.fragment(grammar, derivation)
    @assert frag.state.size !== 0 "Residue $derivation was not found in current grammar. Check if using correct grammar, including NCAAs."
    
    # Remove old sidechain
    for atom in reverse(sidechain)
        ProtoSyn.pop_atom!(pose, atom; keep_downstream_position = false)
        # Note: pop_atom! already re-indexes the pose
    end

    # Insert new sidechain
    frag_sidechain = (!an"^CA$|^N$|^C$|^H$|^O$"r)(frag, gather = true)
    # println("New fragment sidechain: $frag_sidechain")
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
    _ϕ = ProtoSyn.getdihedral(frag.state, frag.graph[1]["C"]) # in template
    ϕ  = ProtoSyn.unit_circle(_ϕ)
    for (index, atom) in enumerate(frag_sidechain) # in template
        parent_is_CA = false
        bonded_to_CA = false
        if atom.parent.name == "CA"
            parent_is_CA     = true
            _atom_dihedral   = ProtoSyn.getdihedral(frag.state, atom)
            atom_dihedral    = ProtoSyn.unit_circle(_atom_dihedral)
            objective_change = atom_dihedral - ϕ
            push!(objective_changes, objective_change)
            ProtoSyn.unbond!(frag, atom, atom.parent)
        elseif "CA" in [a.name for a in atom.bonds]
            bonded_to_CA = true
            ProtoSyn.unbond!(frag, atom, atom.bonds[findfirst(a -> a.name === "CA", atom.bonds)])
        end

        # Insert into the graph
        # Note: insert! already sets the residue.itemsbyname
        insert!(residue, poseCA_index + index, atom)

        # Since now atom is already in the graph, we can bond and add parents
        if parent_is_CA
            ProtoSyn.bond(atom, poseCA)
            setparent!(atom, poseCA)
        elseif bonded_to_CA
            ProtoSyn.bond(atom, poseCA)
        end

        insert!(pose.state, poseCA.index + index, State(T, [frag.state[atom.index]]))
    end

    # Deal with non-contiguous sidechains
    # _start = frag_sidechain[1].index
    # _end   = frag_sidechain[end].index
    # insert!(pose.state, poseCA.index + 1, splice!(frag.state, _start:_end))
    
    reindex(pose.graph)
    reindex(pose.state)

    # Fix CA children positions
    pose_sidechain = (!an"^CA$|^N$|^C$|^H$|^O$"r)(residue, gather = true)
    Δϕ             = pose.state[residue["CA"]].Δϕ
    index          = 1
    _ϕ = ProtoSyn.getdihedral(pose.state, Peptides.phi(residue))
    ϕ  = ProtoSyn.unit_circle(_ϕ)
    for child in residue["CA"].children
        if child in pose_sidechain
            objective = ϕ + objective_changes[index]
            pose.state[child].ϕ = ProtoSyn.unit_circle(objective - Δϕ)
            index += 1
        end
    end
    
    residue.name = ProtoSyn.ResidueName(ProtoSyn.one_2_three[derivation[1][1]])
    ProtoSyn.request_i2c!(pose.state)

    return pose
end

mutate!(pose::Pose{Topology}, sele::AbstractSelection, grammar::LGrammar, derivation; ignore_existing_sidechain::Bool = false) = begin
    residue = sele(pose, gather = true)
    @assert length(residue) == 1 "Can only mutate one residue at a time."
    return mutate!(pose, residue[1], grammar, derivation, ignore_existing_sidechain = ignore_existing_sidechain)
end


"""
    force_mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation, op = "α")

Mutate the given [`Pose`](@ref) `pose` at [`Residue`](@ref) `residue`, changing
it's aminoacid to be `derivation`, as given by the template at `grammar`. This
function changes the whole residue (backbone included). By default, the user
should use [`mutate!`](@ref) instead of this function, except for uncommon
aminoacids.

See also
[`mutate!`](@ref)

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
    phi   = ProtoSyn.getdihedral(pose.state, Peptides.phi(residue))
    
    psi_atom = Peptides.psi(residue)
    if psi_atom !== nothing
        psi = ProtoSyn.getdihedral(pose.state, psi_atom)
    else
        psi = 0.0
    end
    
    omega_atom = Peptides.omega(residue)
    if omega_atom !== nothing
        omega = ProtoSyn.getdihedral(pose.state, omega_atom)
    else
        omega = 0.0
    end

    T     = eltype(pose.state)
    ss    = ProtoSyn.Peptides.SecondaryStructureTemplate{T}(phi, psi, omega)
    old_r_pos = findfirst(residue, residue.container.items)
    Peptides.insert_fragment!(pose, residue, grammar, derivation,
        ss = nothing, op = op)
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
    if downstream_r_pos <= length(ProtoSyn.sequence(pose))
        ProtoSyn.popparent!(new_residue.container.items[downstream_r_pos])
        ProtoSyn.popparent!(new_residue.container.items[downstream_r_pos]["N"])
        
        # 4) Use operator 'op' from 'grammar' to set parenthoods, bonds and correct
        # position, based on the peptidic bond
        grammar.operators[op](new_residue, pose, residue_index = downstream_r_pos)
    end

    # 5) Re-set ascedents
    ProtoSyn.reindex(pose.graph)
    ProtoSyn.request_i2c!(pose.state, all = true)

    return pose
end

# * ----------------------------------------------------------------------------
# * Remove and add sidechains

"""
    remove_sidechains!(pose::Pose{Topology}, res_lib::LGrammar, Opt{AbstractSelection} = nothing)

Removes the sidechain atoms of the given [`Pose`](@ref) `pose`. If an
`AbstractSelection` `selection` is provided, only the sidechain atoms belonging
to the [`Residue`](@ref) instances of that selection are considered for possible
removal. Essentially, the selected [`Residue`](@ref) instances are mutated to
Glycine, based on the provided residue library `res_lib`, without changing the
peptide [`sequence`](@ref ProtoSyn.sequence). Therefore, the original [`sequence`](@ref ProtoSyn.sequence) can be
recovered using the [`add_sidechains!`](@ref) method and energy components such
as [`neighbour_vector`](@ref ProtoSyn.Peptides.Calculators.Caterpillar.neighbour_vector)
can perform correctly.

!!! ukw "Note:"
    Proline residues are ignored.

!!! ukw "Note:"
    Caping [`Atom`](@ref) instances (N-terminal -NH3 and C-terminal -CO2) are
    also removed by this function. See [`cap!`](@ref) to recover from this.

# See also
[`force_remove_sidechains!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Peptides.remove_sidechains!(pose, res_lib)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 147
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)

julia> ProtoSyn.Peptides.remove_sidechains!(pose, res_lib, rn"ALA")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 147
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function remove_sidechains!(pose::Pose{Topology}, res_lib::LGrammar, selection::Opt{AbstractSelection} = nothing)
    _selection = !rn"PRO"
    if selection !== nothing
        _selection = ProtoSyn.promote(_selection & selection, ProtoSyn.Residue)
    end
    selected_residues = _selection(pose, gather = true)
    for residue in selected_residues
        saved_name = residue.name
        Peptides.mutate!(pose, residue, res_lib, seq"G")
        residue.name = saved_name
    end

    sync!(pose)
    return pose
end


"""
    force_remove_sidechains!(pose::Pose{Topology}, selection::Opt{AbstractSelection} = nothing)

Removes the sidechain [`Atom`](@ref) instances of the given [`Pose`](@ref)
`pose`. If an `AbstractSelection` `selection` is provided, only the sidechain
[`Atom`](@ref) instances belonging to that selection are considered for possible
removal. This function completly removes any atom other than backbone atoms (may
break Cα coordination, consider using the [`remove_sidechains!`](@ref) method
instead).

!!! ukw "Note:"
    Proline residues are ignored.

!!! ukw "Note:"
    Caping [`Atom`](@ref) instances (N-terminal -NH3 and C-terminal -CO2) are
    also removed by this function. See [`cap!`](@ref) to recover from this.

# See also
[`remove_sidechains!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Peptides.force_remove_sidechains!(pose)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 105
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)

julia> ProtoSyn.Peptides.force_remove_sidechains!(pose, rn"ALA")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 105
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function force_remove_sidechains!(pose::Pose{Topology}, selection::Opt{AbstractSelection} = nothing)
    _selection = !(an"^CA$|^N$|^C$|^H$|^O$"r | rn"PRO")
    if selection !== nothing
        _selection = _selection & selection
    end
    sidechain = _selection(pose, gather = true)
    for atom in reverse(sidechain) # Note the reverse loop
        ProtoSyn.pop_atom!(pose, atom)
    end

    return pose
end


"""
    add_sidechains!(pose::Pose{Topology}, grammar::LGrammar, selection::Opt{AbstractSelection} = nothing)

Add the sidechain [`Atom`](@ref) instances to the given [`Pose`](@ref) `pose`,
based on the templates of the provided `grammar`. If an `AbstractSelection`
`selection` is given, only the residues of that selection (promoted to `Residue`
instances, using the default aggregator function) are considered for sidechain
addition. The addition is performed using the [`mutate!`](@ref) function, and
follows the current [`Pose`](@ref) [`sequence`](@ref ProtoSyn.sequence).

# Examples
```jldoctest
julia> ProtoSyn.Peptides.add_sidechains!(pose, res_lib)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)

julia> ProtoSyn.Peptides.add_sidechains!(pose, res_lib, rn"ALA")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function add_sidechains!(pose::Pose{Topology}, grammar::LGrammar, selection::Opt{AbstractSelection} = nothing)
    if selection !== nothing
        residues = selection(pose, gather = true)
    else
        residues = collect(eachresidue(pose.graph))
    end
    for residue in residues
        derivation = [string(ProtoSyn.three_2_one[residue.name])]
        derivation == ["P"] && continue
        Peptides.mutate!(pose, residue, grammar, derivation,
            ignore_existing_sidechain = true)
    end

    return pose
end

# * ----------------------------------------------------------------------------
# * Cap and uncap

"""
    is_N_terminal(res::Residue)

Return `true` if the provided [`Residue`](@ref) `residue` is a child of the
residue's container [`root`](@ref).

# Examples
```jldoctest
julia> ProtoSyn.Peptides.is_N_terminal(pose.graph[1][1])
true
```
"""
function is_N_terminal(res::Residue)
    return ProtoSyn.root(res.container.container).container == res.parent
end


"""
    is_C_terminal(res::Residue)

Return `true` if the provided [`Residue`](@ref) `residue` has no children.

# Examples
```jldoctest
julia> ProtoSyn.Peptides.is_C_terminal(pose.graph[1][end])
true
```
"""
function is_C_terminal(res::Residue)
    return length(res.children) == 0
end


"""
    identify_c_terminal(seg::Segment; supress_warn::Bool = false)

Attempts to identify the C terminal of a given [`Segment`](@ref) `seg`, using
the following criteria:
 + Bonding patter: the C terminal follows the bonding pattern C-C-N-C, there the C-terminal is the first C atom of the pattern.
 + Bond number: the C terminal must not exceed 3 bonds.
 + By atom name: the C terminal must be bonded to an atom named CA
As such, this method does not use parenthood relationships to identify the C
terminal. If the `supress_warn` is set to `true` (`false`, by default), any
generated warnings are ignored.

# Examples
```
julia> ProtoSyn.Peptides.identify_c_terminal(pose.graph[1])
Atom{/2a3d:40139/A:1/ASN:73/C:1138}
```
"""
function identify_c_terminal(seg::Segment; supress_warn::Bool = false)

    # Criterium 1. Identify atom by the bonding pattern
    #  - The C terminal is connected to a C, which is connected to a N, which is
    # connected to a C, in the last residue of the segment (based on the residue
    # numbering only)
    last_r = seg[end]
    Cs = ProtoSyn.identify_atom_by_bonding_pattern(last_r, ["C", "C", "N", "C"])
    if isa(Cs, Atom)
        C = Cs
    else
        # Criterium 2. Identify atom by the number of bonds
        #  - The C terminal cannot have more than 3 bonded atoms.
        C = [c for c in Cs if (length(c.bonds) <= 3)]
        if length(C) > 1
            # Criterium 3. Identify atom by a specific bond
            #  - The C terminal should be connected to an atom called "CA"
            C = [c for c in C if c.container["CA"] in c.bonds]
        end
        if length(C) > 1
            !supress_warn && @warn "Multiple candidates for C terminal identification found. Assuming $(C[1]) as C terminal. Check if this is the desired behaviour".
            C = C[1]
        elseif length(C) === 1
            C = C[1]
        elseif length(C) === 0
            !supress_warn && @warn "Terminal C was not found on $seg"
            return Atom[]
        end
    end

    return C
end


"""
    uncap!(pose::Pose, selection::Opt{AbstractSelection} = nothing)

Remove all bonded [`Atom`](@ref) instances to the N- and C-terminal (except Cα)
of the provided [`Pose`](@ref) `pose`. If an `AbstractSelections` `selection` is
provided, search for terminal residues only in the given `selection`. A terminal
is identified based on the following criteria (see [`is_N_terminal`](@ref) and
[`is_C_terminal`](@ref) methods):
- Is a child of the [`Pose`](@ref) `pose` [`root`](@ref);
- Has no children;
Return the modified (in-place) [`Pose`](@ref) `pose`. Does not [`sync!`](@ref).

# Examples
```jldoctest
julia> ProtoSyn.Peptides.uncap!(pose)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 341
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function uncap!(pose::Pose, selection::Opt{AbstractSelection} = nothing; skip_N_terminal::Bool = false, skip_C_terminal::Bool = false)
    # * This function assumes that the N- and C- terminals are correctly named
    # * (N-terminal residue has an atom named "N" and C-terminal has an atom
    # * named "C", both with a bond to a "CA" atom).

    _selection = DownstreamTerminalSelection{Residue}() | UpstreamTerminalSelection{Residue}()
    if selection !== nothing
        _selection = _selection & selection
    end

    residues = _selection(pose, gather = true)
    if length(residues) == 0
        @info "The provided selection/pose doesn't seem to have terminal residues."
        return pose
    end

    for residue in residues
        if is_C_terminal(residue) && !skip_C_terminal
            terminal = residue["C"]
            @assert terminal !== nothing "Can't find atom \"C\" in terminal residue $residue"

            for bond in reverse(terminal.bonds) # Note the reverse loop
                if !(bond.name in ["CA"])
                    ProtoSyn.pop_atom!(pose, bond)
                end
            end
        end

        if is_N_terminal(residue) && !skip_N_terminal # * Note: A residue can be both N- and C- terminal
            terminal = residue["N"]

            for bond in reverse(terminal.bonds) # Note the reverse loop
                if !(bond.name in ["CA"])
                    ProtoSyn.pop_atom!(pose, bond)
                end
            end
        end
    end

    return pose
end


"""
    cap!(pose::Pose, selection::Opt{AbstractSelection} = nothing)

Add template terminals to the N- and C- terminals of the given [`Pose`](@ref)
`pose`. If an `AbstractSelection` `selection` is provided, search for terminal
residues only in the given `selection`. A terminal is identified based on the
following criteria (see [`is_N_terminal`](@ref) and [`is_C_terminal`](@ref)
methods):
- Is a child of the [`Pose`](@ref) `pose` [`root`](@ref);
- Has no children;
Return the modified (in-place) [`Pose`](@ref) `pose`. Performs [`sync!`](@ref)
operation.

# Examples
```jldoctest
julia> ProtoSyn.Peptides.cap!(pose)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 346
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function cap!(pose::Pose, selection::Opt{AbstractSelection} = nothing; skip_N_terminal::Bool = false, skip_C_terminal::Bool = false)
    # * This function assumes that the N- and C- terminals are correctly ordered
    # * (N -> CA -> C, in both terminal backbones).

    # * Default terminal fragments
    T               = eltype(pose.state)
    if !skip_N_terminal
        n_term_filename = Peptides.resource_dir * "/pdb/nterminal.pdb"
        n_term          = ProtoSyn.ProtoSyn.fragment(Peptides.load(T, n_term_filename))
        n_term_atoms    = ["H3", "H2", "H1"] # Reverse order than of desired
    end

    if !skip_C_terminal
        c_term_filename = Peptides.resource_dir * "/pdb/cterminal.pdb"
        c_term          = ProtoSyn.ProtoSyn.fragment(Peptides.load(T, c_term_filename))
        c_term_atoms    = ["OXT", "O"] # Reverse order than of desired
    end

    # * Remove any existing caps
    uncap!(pose, selection,
        skip_N_terminal = skip_N_terminal, skip_C_terminal = skip_C_terminal)
    sync!(pose)

    # * Search for terminals to cap
    _selection = DownstreamTerminalSelection{Residue}() | UpstreamTerminalSelection{Residue}()
    if selection !== nothing
        _selection = _selection & selection
    end
    residues = _selection(pose, gather = true)
    if length(residues) === 0
        @warn "The provided selection/pose doesn't seem to have terminal residues. No capping process will be done. Check if this is the expected behaviour!"
        return pose
    end

    # * Cap terminals
    for residue in residues
        if is_N_terminal(residue) && !skip_N_terminal
            @info "Capping $residue as N-terminal"
            terminal = residue["N"]
            # * Align terminal fragment with existing structure
            mobile_selection = an"^N$|^CA$"r
            target_selection = SerialSelection{Residue}(residue.id, :id) & an"^N$|^CA$"r
            mobile = ProtoSyn.align!(n_term, pose, mobile_selection, target_selection)

            # * Transfer ownership of terminal atoms from fragment to Pose
            # * Placement is based on cartesian coordinates agaisnt the aligned
            # * N terminal template (since internal coordinates depend on the
            # * orientation of the pose in relationship to the root).
            for atom_name in n_term_atoms
                atom        = mobile.graph[1][atom_name]
                atomstate   = mobile.state[atom]
                atom.parent = nothing
                atom.bonds  = Vector{Atom}()
                ProtoSyn.insert_atom_as_children!(pose, terminal, atom, atomstate)
            end

            ProtoSyn.request_c2i!(pose.state)
            sync!(pose)
        end

        if is_C_terminal(residue) && !skip_C_terminal # * Note: A residue can be both N- and C- terminal
            @info "Capping $residue as C-terminal"
            terminal = residue["C"]

            # * Transfer ownership of terminal atoms from fragment to Pose
            # * Placement is based on internal coordinates
            for atom_name in c_term_atoms
                atom        = c_term.graph[1][atom_name]
                atomstate   = c_term.state[atom]
                atom.parent = nothing
                atom.bonds  = Vector{Atom}()
                ProtoSyn.insert_atom_as_children!(pose, terminal, atom, atomstate)
            end

            ProtoSyn.request_i2c!(pose.state)
            sync!(pose)
        end
    end

    return pose
end

include("pose-diagnose.jl")