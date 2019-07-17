@doc raw"""
    apply_ss!(state::State, dihedrals::Vector{Dihedral}, ss::String)

Apply predefined angles to all `PHI` and `PSI` dihedrals defined in `dihedrals`, based on the provided `ss`, changing the State.xyz
to apply the secondary structure. The applied angles (in degrees) are the following:

Beta sheet ("E"):
PHI = -139.0 | PSI = 135.0
Alpha helix ("H"):
PHI = -57.0  | PSI = -47.0

If `override` (Default: `false`) is set to `true`, any previous information saved in `metadata.ss` and `metadata.blocks` is recompiled.

# Examples
```julia-repl
julia> Common.apply_ss!(state, dihedrals, "CCCHHHHCCEEEECCC")
```
See also: [`compile_ss!`](@ref) [`compile_blocks!`](@ref)
"""
function apply_ss!(state::State, metadata::Metadata, ss::String, override::Bool = false)

    # Save secondary structure as metadata
    # if length(metadata.ss) == 0 || override
    #     metadata.ss     = compile_ss!(metadata.dihedrals, ss)
    # end
    if length(metadata.blocks) == 0 || override
        metadata.blocks = compile_blocks!(metadata.residues, ss)
    end

    index::Int64 = 1
    for dihedral in metadata.dihedrals
        #Verify input so that only the phi and psi dihedrals are iterated over, and COILS remain "free"
        if !(dihedral.dtype >= DIHEDRAL.omega || dihedral.residue.ss == Common.SS.COIL)
            rotate_dihedral_to!(state.xyz, dihedral, Common.ss2bbd[dihedral.residue.ss][dihedral.dtype])
        end
    end
    printstyled("(SETUP) ▲ Applied secondary structure angles to $(length(metadata.blocks)) blocks\n", color = 9)
end


@doc raw"""
    apply_backbone_from_file!(target_state::State, target_dihedrals::Vector{Dihedral}, template_file::String)

Read the input file `template_file`, and apply both dihedral and regular angles to the backbone of `target_state`.
`target_dihedrals` should contain backbone dihedrals only (phi, psi, omega).
Note: `template_file` needs to have CONECT records.

# Examples
```julia-repl
julia> Common.apply_backbone_dihedrals_from_file(state, bb_dihedrals, "native_conf.pdb")
```
"""
function apply_backbone_from_file!(target_state::State, target_dihedrals::Vector{Dihedral}, template_file::String)
    Common.apply_backbone_angles_from_file!(target_state, target_dihedrals, template_file)
    Common.apply_backbone_dihedrals_from_file!(target_state, target_dihedrals, template_file)
end

@doc raw"""
    apply_backbone_dihedrals_from_file!(target_state::State, target_dihedrals::Vector{Dihedral}, template_file::String)

Read the input file `template_file`, extract the dihedrals angles and apply to the backbone dihedrals `target_dihedrals` of `target_state`.
`target_dihedrals` should contain backbone dihedrals only (phi, psi, omega).
Note: `template_file` needs to have CONECT records.

# Examples
```julia-repl
julia> Common.apply_backbone_dihedrals_from_file(state, bb_dihedrals, "native_conf.pdb")
```
"""
function apply_backbone_dihedrals_from_file!(target_state::State, target_dihedrals::Vector{Dihedral}, template_file::String)

    template_state, template_metadata = Common.load_from_pdb(template_file)
    template_dihedrals = filter(x -> x.dtype <= Common.DIHEDRAL.omega, template_metadata.dihedrals)
    for (cd, template_dihedral) in enumerate(template_dihedrals)
        target_angle = Aux.calc_dihedral(target_state.xyz[target_dihedrals[cd].a1, :], target_state.xyz[target_dihedrals[cd].a2, :], target_state.xyz[target_dihedrals[cd].a3, :], target_state.xyz[target_dihedrals[cd].a4, :])
        template_angle = Aux.calc_dihedral(template_state.xyz[template_dihedral.a1, :], template_state.xyz[template_dihedral.a2, :], template_state.xyz[template_dihedral.a3, :], template_state.xyz[template_dihedral.a4, :])
        displacement = template_angle - target_angle
        rotate_dihedral!(target_state.xyz, target_dihedrals[cd], displacement)
    end
    printstyled("(SETUP) ▲ Applied backbone dihedrals to $(length(template_dihedrals)) dihedrals\n", color = 9)
end


@doc raw"""
    apply_backbone_angles_from_file!(target_state::State, target_dihedrals::Vector{Dihedral}, template_file::String)

Read the input file `template_file`, extract the backbone angles and apply to the backbone atoms of `target_state`.
`target_dihedrals` should contain backbone dihedrals only (phi, psi, omega).
Note: `template_file` needs to have CONECT records.

# Examples
```julia-repl
julia> Common.apply_backbone_angles_from_file(state, bb_dihedrals, "native_conf.pdb")
```
"""
function apply_backbone_angles_from_file!(target_state::State, target_dihedrals::Vector{Dihedral}, template_file::String)

    template_state, template_metadata = Common.load_from_pdb(template_file)
    template_dihedrals = filter(x -> x.dtype <= Common.DIHEDRAL.omega, template_metadata.dihedrals)
    for (cd, template_dihedral) in enumerate(template_dihedrals)
        target_angle = Aux.calc_angle(target_state.xyz[target_dihedrals[cd].a1, :], target_state.xyz[target_dihedrals[cd].a2, :], target_state.xyz[target_dihedrals[cd].a3, :])
        template_angle = Aux.calc_angle(template_state.xyz[template_dihedral.a1, :], template_state.xyz[template_dihedral.a2, :], template_state.xyz[template_dihedral.a3, :])
        displacement = template_angle - target_angle

        v21   = target_state.xyz[target_dihedrals[cd].a1, :] - target_state.xyz[target_dihedrals[cd].a2, :]
        v23   = target_state.xyz[target_dihedrals[cd].a3, :] - target_state.xyz[target_dihedrals[cd].a2, :]
        pivot = target_state.xyz[target_dihedrals[cd].a2, :]
        rmat  = Aux.rotation_matrix_from_axis_angle(cross(v21, v23), displacement)
        movable = collect(target_dihedrals[cd].a2:size(target_state.xyz, 1))
        target_state.xyz[movable, :] = ((rmat * (target_state.xyz[movable, :] .- pivot')') .+ pivot)'
    end
    printstyled("(SETUP) ▲ Applied backbone angles to $(length(template_dihedrals)) angles\n", color = 9)
end



@doc raw"""
    stretch_conformation!(state::State, dihedrals::Vector{Dihedral})

Iterate over the `dihedrals` list and rotate all PHI and PSI to 180º.

# Examples
```julia-repl
julia> Common.stretch_conformation!(state, dihedrals)
```
"""
function stretch_conformation!(state::State, dihedrals::Vector{Dihedral})

    angles = map(x -> deg2rad(x), [111.099, 119.651, 117.817])
    index = 1
    for dihedral in dihedrals
        rotate_dihedral_to!(state.xyz, dihedral, 3.1416)

        target_angle = angles[index]
        current_angle = Aux.calc_angle(state.xyz[dihedral.a1, :], state.xyz[dihedral.a2, :], state.xyz[dihedral.a3, :])
        displacement = target_angle - current_angle

        v21   = state.xyz[dihedral.a1, :] - state.xyz[dihedral.a2, :]
        v23   = state.xyz[dihedral.a3, :] - state.xyz[dihedral.a2, :]
        pivot = state.xyz[dihedral.a2, :]
        rmat  = Aux.rotation_matrix_from_axis_angle(cross(v21, v23), displacement)
        movable = collect(dihedral.a2:size(state.xyz, 1))
        state.xyz[movable, :] = ((rmat * (state.xyz[movable, :] .- pivot')') .+ pivot)'

        index += 1
        if index > length(angles)
            index = 1
        end
    end
end

#-----------------------------------#
# L A T T I C E   G E N E R A T O R #
#-----------------------------------#

module LATTICE
    @enum TYPE begin
        primitive       = 1
        body_centered   = 2
        face_centered   = 3
    end
end

function generate_template(hf::Vector{Float64}, type::LATTICE.TYPE)::Array{Float64, 2}
    atoms = [0.0 0.0 0.0]
    if type == LATTICE.body_centered
        atoms = vcat(atoms, [hf[1] hf[2] hf[3]])
    elseif type == LATTICE.face_centered
        atoms = vcat(atoms, [hf[1] hf[2] 0.0; hf[1] 0.0 hf[3]; 0.0 hf[2] hf[3]])
    end
    return atoms
end

function calculate_n_atoms(rep::Vector{Int64}, type::LATTICE.TYPE,
    closed::Bool = true)::Int64

    if type == LATTICE.primitive
        n_atoms = (rep[1] + closed) * (rep[2] + closed) * (rep[3] + closed)
    elseif type == LATTICE.body_centered
        n_atoms = (rep[1] + closed) * (rep[2] + closed) * (rep[3] + closed)
        n_atoms += rep[1]*rep[2]*rep[3]
    elseif type == LATTICE.face_centered
        n_atoms = (rep[1] + closed) * (rep[2] + closed) * (rep[3] + closed)
        n_atoms += rep[1]*rep[2]*rep[3] * 3
        if closed
            n_atoms += rep[1]*rep[2] + rep[2] * rep[3] + rep[3] * rep[1]
        end
    end
    return n_atoms
end

function generate_cubic_lattice(size::Float64, rep::Int64,
    type::LATTICE.TYPE, closed::Bool = true)::State
    return generate_lattice([size, size, size], [rep, rep, rep], type, closed)
end

function generate_lattice(side_length::Vector{Float64}, rep::Vector{Int64},
    type::LATTICE.TYPE, closed::Bool = true)::State
    # Side length in nm

    # Generate the template
    f = (side_length./rep)
    template = generate_template(f./2, type)

    # Create empty list of atoms
    n_atoms  = calculate_n_atoms(rep, type, closed)
    state    = State(n_atoms)

    # Fill State with repetitions of template
    c = 1
    for i in 0:rep[1]
        for j in 0:rep[2]
            for k in 0:rep[3]
                for index in 1:size(template, 1)
                    new_pos = template[index, :] .+ [i*f[1], j*f[2], k*f[3]]
                    if closed && all(new_pos .<= side_length) || all(new_pos .< side_length)
                        state.xyz[c, :] = new_pos'
                        c +=1
                    end
                end
            end
        end
    end
    return state
end


