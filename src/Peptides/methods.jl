export grammar

"""
    grammar([::Type{T}];[verbose::Bool = true]) where {T <: AbstractFloat}

Build a `LGrammar` for peptides, taking as variables the fragments in the
default resource directory. If the option type T is not provided, the default
ProtoSyn float value will be used. The returned LGrammar is required for
building peptides from fragments. If `verbose` is set to `true` (is, by
default), print the loading status.

# Examples
```jldoctest
julia> pose = ProtoSyn.Peptides.build(res_lib, seq"AAGASTASSE")
┌ Warning: Residue Residue{/UNK:1/UNK:1/GLU:10} has no psi angle
└ @ ProtoSyn.Peptides.Dihedral ~/project_c/ProtoSyn.jl/src/Peptides/constants.jl:67
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 109
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function grammar(::Type{T}; verbose::Bool = true) where {T <: AbstractFloat}
    filename = joinpath(Peptides.resource_dir, "grammars.yml")
    ProtoSyn.load_grammar_from_file(T, filename, "peptide", verbose = verbose)
end

grammar(;verbose::Bool = true) = grammar(ProtoSyn.Units.defaultFloat; verbose = verbose)


"""
    build(grammar::LGrammar{T}, derivation, ss::NTuple{3,Number} = SecondaryStructure[:linear]) where {T <: AbstractFloat}

Build a `Pose{Topology}` using the given `derivation` sequence on the provided
`grammar` instructions. If an `ss` is provided, automatically apply it to the
built pose (linear secondary structure, by default).
!!! note
    This function is an overload of `ProtoSyn.build`.

# See also
`setss!`

# Examples
```jldoctest
julia> pose = ProtoSyn.Peptides.build(res_lib, seq"QQQ")
┌ Warning: Residue Residue{/UNK:1/UNK:1/GLN:3} has no psi angle
└ @ ProtoSyn.Peptides.Dihedral ~/project_c/ProtoSyn.jl/src/Peptides/constants.jl:67
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 51
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function build(grammar::LGrammar{T}, derivation, ss::NTuple{3,Number} = SecondaryStructure[:linear]) where {T <: AbstractFloat}

    pose = ProtoSyn.build(grammar, derivation)
    ProtoSyn.Peptides.setss!(pose, ss)
    sync!(pose)
    pose
end


"""
    append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ss::NTuple{3,Number} = SecondaryStructure[:linear], op = "α")

Based on the provided `grammar`, add the residue sequence from `derivation` to
the given `pose`, appending it AFTER the given `residue`. This residue and the
new fragment will be connected using operation `op` ("α" by default). Set the
secondary structure of the added residue to match the given `ss` (linear, by
default). Return the altered `pose`.

# Examples
```jldoctest
julia> ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][1], res_lib, seq"A")
┌ Warning: Residue Residue{/UNK:1/UNK:1/ALA:2} has no psi angle
└ @ ProtoSyn.Peptides.Dihedral ~/project_c/ProtoSyn.jl/src/Peptides/constants.jl:67
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 353
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function append_fragment!(pose::Pose{Topology}, residue::Residue,
    grammar::LGrammar, derivation;
    ss::NTuple{3,Number} = SecondaryStructure[:linear], op = "α")

    ProtoSyn.append_fragment!(pose, residue, grammar, derivation; op = op)
    residues = residue.container.items[(residue.index + 1):(residue.index + length(derivation))]
    setss!(pose, ss, residues)
    return pose
end

"""
    Careful: Rotamer. 
"""
function append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, frag::Pose{Segment}; ss::Opt{NTuple{3,Number}} = nothing, op = "α")
    
    residue_selection = SerialSelection{Residue}(residue.id, :id)
    rotamer = Rotamers.get_rotamer(pose, residue)
    ProtoSyn.append_fragment!(pose, residue, grammar, frag, op = op)
    
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
julia> ProtoSyn.Peptides.insert_residues!(pose, pose.graph[1][2], res_lib, seq"A")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 353
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function insert_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation;
    ss::NTuple{3, Number} = SecondaryStructure[:linear], op = "α")

    connected_to_origin = residue.parent == ProtoSyn.root(pose.graph).container
    if connected_to_origin
        state = pose.state
        N = state[residue["N"]] # This is the N atom state
        i = residue.index       # This is the residue index
        (b, θ, ϕ) = (N.b, N.θ, N.ϕ)
    else
        ProtoSyn.unbond!(pose, residue.container[residue.index - 1]["C"], residue["N"])
    end
    
    ProtoSyn.insert_fragment!(pose, residue, grammar, derivation; op = op,
        connect_upstream = !connected_to_origin)

    if connected_to_origin
        N = state[residue.container[i]["N"]]
        (N.b, N.θ, N.ϕ) = (b, θ, ϕ)
        ProtoSyn.request_i2c!(state; all = true)
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
julia> ProtoSyn.Peptides.force_mutate!(pose, pose.graph[1][3], res_lib, seq"K")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 343
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
    old_r_pos = findfirst(residue, residue.container.items)
    Peptides.insert_residues!(pose, residue, grammar, derivation,
        ss = (phi, psi, omega), op = op)
    new_residue = residue.container.items[old_r_pos]

    # 2) Remove old residue (now in position R + 1)
    old_r_pos += 1
    N = residue.container.items[old_r_pos]["N"]
    Peptides.pop_residue!(pose, residue.container.items[old_r_pos])
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
        ProtoSyn.Peptides.pop_atom!(pose, atom)
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
    pop_residue!(pose::Pose{Topology}, residue::Residue)

Delete `residue` from the given `pose`. Removes parenthood (at residue and atom
level), removes bonds, reindexes ascedents and sets position and parenthood to
origin of downstream residues.

# Examples
```
julia> ProtoSyn.Peptides.pop_residue!(pose, pose.graph[1][3])
Pose{Residue}(Residue{/SER:16610}, State{Float64}:
 Size: 11
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function pop_residue!(pose::Pose{Topology}, residue::Residue)

    # 1) Since we are calling Peptides, we know that we should unbond this
    # residue and the next (children).
    # * Note: Peptides.unbond!sets the position of the child residue, while
    # * connecting it to the root
    for child in residue.children
        Peptides.unbond!(pose, residue, child) # Order is important
    end

    # 2) Now we can safelly pop the residue, while maintaining the positions of
    # downstream residues
    ProtoSyn.pop_residue!(pose, residue)
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
    Caping atoms (N-terminal -NH3 and C-terminal -CO2) are also removed by this function. See [`cap!`](@ref) to recover from this.

# See also
`force_remove_sidechains!`

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
    Caping atoms (N-terminal -NH3 and C-terminal -CO2) are also removed by this function. See [`cap!`](@ref) to recover from this.

# See also
`remove_sidechains!`

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
        ProtoSyn.Peptides.pop_atom!(pose, atom)
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
        derivation = [string(Peptides.three_2_one[residue.name])]
        derivation == ["P"] && continue
        Peptides.mutate!(pose, residue, grammar, derivation,
            ignore_existing_sidechain = true)
    end

    return pose
end


"""
    is_N_terminal(res::Residue)

Return `true` if the provided `residue` is a child of the residue's container.

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

Return `true` if the provided `residue` has no children.

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
    uncap!(pose::Pose, selection::Opt{AbstractSelection} = nothing)

Remove all bonded atoms to the N- and C-terminal (except Cα) of the provided
`pose`. If a `selection` is provided, search for terminal residues only in the
given residue. A terminal is identified based on the following criteria:
- Is a child of the `pose` origin;
- Has no children;
Return the modified (in-place) `pose`. Does not sync!

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
function uncap!(pose::Pose, selection::Opt{AbstractSelection} = nothing)
    # * This function assumes that the N- and C- terminals are correctly named
    # * (N-terminal residue has an atom named "N" and C-terminal has an atom
    # * named "C", both with a bond to a "CA" atom).

    _selection = ProtoSyn.TerminalSelection()
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
                    ProtoSyn.Peptides.pop_atom!(pose, bond)
                end
            end
        end

        if is_N_terminal(residue) # * Note: A residue can be both N- and C- terminal
            terminal = residue["N"]

            for bond in reverse(terminal.bonds) # Note the reverse loop
                if !(bond.name in ["CA"])
                    ProtoSyn.Peptides.pop_atom!(pose, bond)
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
julia> ProtoSyn.Peptides.cap!(pose)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function cap!(pose::Pose, selection::Opt{AbstractSelection} = nothing)
    # * This function assumes that the N- and C- terminals are correctly ordered
    # * (N -> CA -> C, in both terminal backbones).

    # * Default terminal fragments
    T               = eltype(pose.state)
    n_term_filename = Peptides.resource_dir * "/pdb/nterminal.pdb"
    n_term          = ProtoSyn.ProtoSyn.fragment(Peptides.load(T, n_term_filename; verbose = false))
    n_term_atoms    = ["H1", "H2", "H3"]

    c_term_filename = Peptides.resource_dir * "/pdb/cterminal.pdb"
    c_term          = ProtoSyn.ProtoSyn.fragment(Peptides.load(T, c_term_filename; verbose = false))
    c_term_atoms    = ["O", "OXT"]

    # * Remove any existing caps
    uncap!(pose, selection)

    # * Search for terminals to cap
    _selection = ProtoSyn.TerminalSelection()
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
                ProtoSyn.insert_atom_as_children!(pose, terminal, atom, atomstate)
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
                ProtoSyn.insert_atom_as_children!(pose, terminal, atom, atomstate)
            end
        end
    end

    return pose
end


"""
    # TODO
"""
function pop_atom!(pose::Pose{Topology}, atom::Atom)::Pose{Atom}
    ProtoSyn.pop_atom!(pose, atom; unbond_f = ProtoSyn.Peptides.unbond!)
end