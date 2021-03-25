export grammar

"""
    grammar([::Type{T}]) where {T <: AbstractFloat}

Build a `LGrammar` for peptides, taking as variables the fragments in the
default resource directory. If the option type T is not provided, the default
ProtoSyn float value will be used. The returned LGrammar is required for
building peptides from fragments.

# Examples
```julia-repl
julia> g = Peptides.grammar();
julia> pose = Peptides.build(grammar, seq"AAGASTASSE")
...
```
"""
function grammar(::Type{T}) where {T <: AbstractFloat}
    filename = joinpath(Peptides.resource_dir, "grammars.yml")
    Builder.fromfile(T, filename, "peptide")
end

grammar() = grammar(ProtoSyn.Units.defaultFloat)


"""
    build(grammar::LGrammar{T}, derivation, ss::NTuple{3,Number} = SecondaryStructure[:linear]) where {T <: AbstractFloat}

Build a `Pose{Topology}` using the given `derivation` sequence on the provided
`grammar` instructions. If an `ss` is provided, automatically apply it to the
built pose (linear secondary structure, by default).
!!! note
    This function is an overload of `Builder.build`.

# See also
`setss!`

# Examples
```julia-repl
julia> pose = Peptides.build(grammar, seq"QQQ");
...
```
"""
function build(grammar::LGrammar{T}, derivation, ss::NTuple{3,Number} = SecondaryStructure[:linear]) where {T <: AbstractFloat}

    pose = Builder.build(grammar, derivation)
    Peptides.setss!(pose, ss)
    sync!(pose)
    pose
end

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
            Dihedral.psi(r) !== nothing && ProtoSyn.setdihedral!(container.state, Dihedral.psi(r),  ψ)
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

function append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, frag::Pose{Segment}; ss::Opt{NTuple{3,Number}} = nothing, op = "α")
    
    residue_selection = SerialSelection{Residue}(residue.id, :id)
    rotamer = Rotamers.get_rotamer(pose, residue)
    Builder.append_fragment!(pose, residue, grammar, frag, op = op)
    
    if ss !== nothing
        residues = residue.container.items[(residue.index + 1):(residue.index + length(frag.graph))]
        setss!(pose, ss, residues)
    end

    derivation = [string(Peptides.three_2_one[residue.name.content])]
    Peptides.force_mutate!(pose, residue, grammar, derivation)
    residue = residue_selection(pose, gather = true)[1]
    Rotamers.apply!(pose.state, rotamer, residue)
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
    N = residue.container.items[old_r_pos]["N"]
    Peptides.pop_residue!(pose, residue.container.items[old_r_pos])
    # i = findall(a -> a == N, ProtoSyn.origin(pose.graph).children)
    # deleteat!(ProtoSyn.origin(pose.graph).children, i)

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
    return pose
end


"""
    mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation)

Mutate the given `pose` at `residue`, changing it's aminoacid to be
`derivation`, as given by the template at `grammar`. This function changes the
sidechain only (± 7x faster than `force_mutate!`). When mutating to Proline,
falls back to `force_mutate!`. If `ignore_existing_sidechain` is set to `true`
(`false` by default), existing sidechains are first removed and then re-added,
regardless of being of the same type.
!!! note
    Sidechains are selected based on the atom name (backbone atoms must be named N, H, CA, C and O, exclusively. Non backbone atoms should have other names, such as H1, H2, etc.)

# Examples
```jldoctest
julia> Peptides.mutate!(pose, pose.graph[1][3], res_lib, seq"K")
```
"""
function mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ignore_existing_sidechain::Bool = false)

    @assert length(derivation) == 1 "Derivation must have length = 1."
    
    sidechain = (!an"^CA$|^N$|^C$|^H$|^O$"r)(residue, gather = true)

    same_aminoacid = string(Peptides.three_2_one[residue.name]) == derivation[1]
    if same_aminoacid && length(sidechain) > 0 && !ignore_existing_sidechain
        # println("No mutation required, residue already has sidechain of the requested type.")
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
    # residue and the next (children).
    # * Note: Peptides.unbond sets the position of the child residue, while
    # * connecting it to the root
    for child in residue.children
        Peptides.unbond(pose, residue, child) # Order is important
    end

    # 2) Now we can safelly pop the residue, while maintaining the positions of
    # downstream residues
    Builder.pop_residue!(pose, residue)
end


"""
    remove_sidechains!(pose::Pose{Topology}, res_lib::selection::LGrammar, Opt{AbstractSelection} = nothing)

Removes the sidechain atoms of the given `pose`. If a `selection` is given, only
the sidechain atoms belonging to the residues of that selection are considered
for possible removal. Essentially, the selected residues are mutated to Glycine,
based on the provided residue library `res_lib`, without changing the peptide
sequence. Therefore, the original sequence can be recovered using
`Peptides.add_sidechains!` and energy components such as `calc_solvation_energy`
can perform correctly.
!!! note
    Proline residues are ignored.
!!! note
    Caping atoms (N-terminal -NH3 and C-terminal -CO2) are also removed by this function. See `cap!` to recover from this.

# See also
`force_remove_sidechains!`

# Examples
```jldoctest
julia> Peptides.remove_sidechains!(pose, res_lib)

julia> Peptides.remove_sidechains!(pose, res_lib, rn"ALA")
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

    return pose
end


"""
    force_remove_sidechains!(pose::Pose{Topology}, res_lib::selection::LGrammar, Opt{AbstractSelection} = nothing)

Removes the sidechain atoms of the given `pose`. If a `selection` is given, only
the sidechain atoms belonging to the residues of that selection are considered
for possible removal. This function completly removes any atom other than
backbone atoms (may break Cα coordination).

!!! note
    Proline residues are ignored.
!!! note
    Caping atoms (N-terminal -NH3 and C-terminal -CO2) are also removed by this function. See `cap!` to recover from this.

# See also
`remove_sidechains!`

# Examples
```jldoctest
julia> Peptides.force_remove_sidechains!(pose)

julia> Peptides.force_remove_sidechains!(pose, rn"ALA")
```
"""
function force_remove_sidechains!(pose::Pose{Topology}, selection::Opt{AbstractSelection} = nothing)
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
        derivation == ["P"] && continue
        Peptides.mutate!(pose, residue, grammar, derivation,
            ignore_existing_sidechain = true)
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


"""
    is_N_terminal(res::Residue)

Return `true` if the provided `residue` is a child of the residue's container.

# Examples
```jldoctest
julia> Peptides.is_N_terminal(pose.graph[1][1])
true
```
"""
function is_N_terminal(res::Residue)
    return ProtoSyn.origin(res.container.container).container == res.parent
end


"""
    is_C_terminal(res::Residue)

Return `true` if the provided `residue` has no children.

# Examples
```jldoctest
julia> Peptides.is_N_terminal(pose.graph[1][end])
true
```
"""
function is_C_terminal(res::Residue)
    return length(res.children) == 0
end


"""
    uncap!(pose::Pose, selection::Opt{AbstractSelection} = nothing)

Remove all bonded atoms to the N- and C-terminal (except Cα) of the provided
`pose`. If a `selection` is provided, search for terminal residues only in the
given residue. A terminal is identified based on the following criteria:
- Is a child of the `pose` origin;
- Has no children;
Return the modified (in-place) `pose`.

# Examples
```jldoctest
julia> Peptides.uncap!(pose)
```
"""
function uncap!(pose::Pose, selection::Opt{AbstractSelection} = nothing)
    # * This function assumes that the N- and C- terminals are correctly named
    # * (N-terminal residue has an atom named "N" and C-terminal has an atom
    # * named "C", both with a bond to a "CA" atom).

    _selection = ProtoSyn.TerminalSelection{Residue}()
    if selection !== nothing
        _selection = _selection & selection
    end

    residues = _selection(pose, gather = true)
    if length(residues) == 0
        @info "The provided selection/pose doesn't seem to have terminal residues."
        return pose
    end

    for residue in residues
        if is_C_terminal(residue)
            terminal = residue["C"]

            for bond in reverse(terminal.bonds) # Note the reverse loop
                if !(bond.name in ["CA"])
                    Builder.pop_atom!(pose, bond)
                end
            end
        end

        if is_N_terminal(residue) # * Note: A residue can be both N- and C- terminal
            terminal = residue["N"]

            for bond in reverse(terminal.bonds) # Note the reverse loop
                if !(bond.name in ["CA"])
                    Builder.pop_atom!(pose, bond)
                end
            end
        end
    end

    return pose
end


"""
    cap!(pose::Pose, selection::Opt{AbstractSelection} = nothing)

Add template terminals to the N- and C- terminals. If a `selection` is provided,
search for terminal residues only in the given residue. A terminal is identified
based on the following criteria:
- Is a child of the `pose` origin;
- Has no children;
Return the modified (in-place) `pose`.

# Examples
```jldoctest
julia> Peptides.cap!(pose)
```
"""
function cap!(pose::Pose, selection::Opt{AbstractSelection} = nothing)
    # * This function assumes that the N- and C- terminals are correctly ordered
    # * (N -> CA -> C, in both terminal backbones).

    # * Default terminal fragments
    T               = eltype(pose.state)
    n_term_filename = Peptides.resource_dir * "/pdb/nterminal.pdb"
    n_term          = ProtoSyn.Builder.fragment(Peptides.load(T, n_term_filename))
    n_term_atoms    = ["H1", "H2", "H3"]

    c_term_filename = Peptides.resource_dir * "/pdb/cterminal.pdb"
    c_term          = ProtoSyn.Builder.fragment(load(T, c_term_filename))
    c_term_atoms    = ["O", "OXT"]

    # * Remove any existing caps
    uncap!(pose, selection)

    # * Search for terminals to cap
    _selection = ProtoSyn.TerminalSelection{Residue}()
    if selection !== nothing
        _selection = _selection & selection
    end
    residues = _selection(pose, gather = true)
    @assert length(residues) > 0 "The provided selection/pose doesn't seem to have terminal residues."

    # * Cap terminals
    for residue in residues
        if is_N_terminal(residue)
            terminal = residue["N"]
            # * Align terminal fragment with existing structure
            mobile_selection = an"^N$|^CA$|^C$"r
            target_selection = SerialSelection{Residue}(residue.id, :id) & an"^N$|^CA$|^C$"r
            mobile = ProtoSyn.align!(n_term, pose, mobile_selection, target_selection)

            # * Transfer ownership of terminal atoms from fragment to Pose
            for atom_name in n_term_atoms
                atom        = mobile.graph[1][atom_name]
                atomstate   = mobile.state[atom]
                atom.parent = nothing
                atom.bonds  = Vector{Atom}()
                Builder.insert_atom_as_children!(pose, terminal, atom, atomstate)
            end
        end
        if is_C_terminal(residue) # * Note: A residue can be both N- and C- terminal
            terminal = residue["C"]
            # * Align terminal fragment with existing structure
            mobile_selection = an"^N$|^CA$|^C$"r
            target_selection = SerialSelection{Residue}(residue.id, :id) & an"^N$|^CA$|^C$"r
            mobile = ProtoSyn.align!(c_term, pose, mobile_selection, target_selection)

            # * Transfer ownership of terminal atoms from fragment to Pose
            for atom_name in c_term_atoms
                atom        = mobile.graph[1][atom_name]
                atomstate   = mobile.state[atom]
                atom.parent = nothing
                atom.bonds  = Vector{Atom}()
                Builder.insert_atom_as_children!(pose, terminal, atom, atomstate)
            end
        end
    end
end


"""
    sequence(container::ProtoSyn.AbstractContainer)::String

Return the sequence of aminoacids (in 1 letter mode) of the given container/pose
as a string.

# Examples
```julia-repl
julia> sequence(pose)
"AAGASTASSE"
```
"""
function sequence(container::ProtoSyn.AbstractContainer)::String

    sequence = ""
    for residue in eachresidue(container)
        try
            sequence *= three_2_one[residue.name]
        catch KeyError
            sequence *= '?'
        end
    end

    return sequence
end

sequence(pose::Pose) = sequence(pose.graph)