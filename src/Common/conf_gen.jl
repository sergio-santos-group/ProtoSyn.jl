@doc raw"""
    apply_ss!(state::State, dihedrals::Vector{Dihedral}, ss::String)

Apply predefined angles to all `PHI` and `PSI` dihedrals defined in `dihedrals`, based on the provided `ss`, changing the State.xyz
to apply the secondary structure. The applied angles (in degrees) are the following:

Beta sheet ("E"):
PHI = -139.0 | PSI = 135.0
Alpha helix ("H"):
PHI = -57.0  | PSI = -47.0

# Examples
```julia-repl
julia> Common.apply_ss!(state, dihedrals, "CCCHHHHCCEEEECCC")
```
See also: [`infer_ss`](@ref)
"""
function apply_ss!(state::State, dihedrals::Vector{Dihedral}, ss::String)

    # Save secondary structure as metadata
    state.metadata.ss = infer_ss(dihedrals, ss)

    # Define target values for dihedral angles
    t_angles = Dict(
        'E' => Dict(phi => deg2rad(-139.0), psi => deg2rad(135.0)),
        'H' => Dict(phi => deg2rad(-57.0),  psi => deg2rad(-47.0)))

    index::Int64 = 0
    for dihedral in dihedrals
        #Verify input so that only the phi and psi dihedrals are iterated over
        dihedral.dtype > omega ? continue : nothing
        dihedral.dtype == psi ? index += 1 : nothing
        index > length(ss) ? error("The length of the secondary stucture string is inferior to the phi and psi count in the molecule.") : nothing
        ss[index] == 'C' ? continue : nothing

        rotate_dihedral_to!(state.xyz, dihedral, t_angles[ss[index]][dihedral.dtype])
    end

end


@doc raw"""
    infer_ss(dihedrals::Vector{Dihedral}, ss::String)::Vector{SecondaryStructureMetadata}
"""
function infer_ss(dihedrals::Vector{Dihedral}, ss::String)::Vector{SecondaryStructureMetadata}

    residues = Vector{Residue}()
    for dihedral in dihedrals
        if !(dihedral.residue in residues)
            push!(residues, dihedral.residue)
        end
    end
    return infer_ss(residues, ss)
end

@doc raw"""
    infer_ss(residues::Dict{Int64, Residue}, ss::String)::Vector{SecondaryStructureMetadata}
"""
function infer_ss(residues::Dict{Int64, Residue}, ss::String)::Vector{SecondaryStructureMetadata}
    infer_ss(collect(values(residues)), ss)
end

@doc raw"""
    infer_ss(residues::Vector{Residue}, ss::String)::Vector{SecondaryStructureMetadata}

Read the provided secondary structure string `ss` and infer the [`SecondaryStructureMetadata`](@ref) [`Metadata`](@ref) from the provided list of residues/dihedrals.

# Examples
```julia-repl
julia> Common.infer_ss(dihedrals, "CCCHHHHCCEEEECCC")
2-element Array{ProtoSyn.Common.SecondaryStructureMetadata,1}:
 SecondaryStructureMetadata(ss_type=HELIX, name=HA, I-4 <-> A-7, conf=1)  
 SecondaryStructureMetadata(ss_type=SHEET, name=BA, A-10 <-> V-13, conf=1)
```
"""
function infer_ss(residues::Vector{Residue}, ss::String)::Vector{SecondaryStructureMetadata}

    conv_type = Dict('H' => SS.HELIX, 'E' => SS.SHEET, 'C' => SS.COIL)
    conv_name = Dict('H' => "HA", 'E' => "BA")

    sec_str = SecondaryStructureMetadata[]
    last_ss::Char = ss[1]
    i_idx::Int64 = 1
    for (index, curr_ss) in enumerate(ss)
        if curr_ss != last_ss
            if last_ss in ['H', 'E']
                push!(sec_str, SecondaryStructureMetadata(conv_type[last_ss], conv_name[last_ss], residues[i_idx].name, i_idx, residues[index - 1].name, index - 1, 1))
            end
            i_idx = index
        end
        residues[index].ss = conv_type[curr_ss] #If commented, will not apply ss to residue.ss
        last_ss = curr_ss
    end
    return sec_str
end


@doc raw"""
    apply_dihedrals_from_file!(state::State, bb_dihedrals::Vector{Dihedral}, file_i::String)

Read the input file `file_i`, extract the dihedrals angles and apply to the backbone dihedrals `bb_dihedrals` of `state`.

# Examples
```julia-repl
julia> Common.apply_initial_conf(state, bb_dihedrals, "native_conf.pdb")
```
"""
function apply_dihedrals_from_file!(state::State, bb_dihedrals::Vector{Dihedral}, file_i::String)

    #Extract backbone coordinates
    backbone = Vector{Vector{Float64}}()
    open(file_i, "r") do f
        for line in eachline(f)
            if startswith(line, "ATOM") && string(strip(line[13:16])) in ["N", "CA", "C"]
                push!(backbone, map(x -> 0.1*parse(Float64, x), [line[31:38], line[39:46], line[47:54]]))
            end
        end
    end

    #Calculate native dihedrals and rotate on state
    count::Int64 = 0
    cd = 1
    for i in 1:(length(backbone) - 3)
        if i == (2 + 3 * count)
            count += 1
            continue
        end
        current_angle = Aux.calc_dih_angle(state.xyz[bb_dihedrals[cd].a1, :], state.xyz[bb_dihedrals[cd].a2, :], state.xyz[bb_dihedrals[cd].a3, :], state.xyz[bb_dihedrals[cd].a4, :])
        displacement = Aux.calc_dih_angle(backbone[i], backbone[i + 1], backbone[i + 2], backbone[i + 3]) - current_angle
        rotate_dihedral!(state.xyz, bb_dihedrals[cd], displacement)
        cd += 1
    end

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

    for dihedral in dihedrals
        dihedral.dtype > omega ? continue : nothing
        rotate_dihedral_to!(state.xyz, dihedral, 3.1416)
    end
end



@doc raw"""
    fix_proline!(state::State, dihedrals::Vector{Dihedral})

( TEMPORARY ) Find all Prolines in the `dihedrals` list and solve the two following issues:
1. Add the first two atoms to the movables list.
2. Rotate the dihedral to be at 0º instead of 180º.

# Examples
```julia-repl
julia> Common.fix_proline(state, dihedrals)
```
"""
function fix_proline!(state::State, dihedrals::Vector{Dihedral})

    for dihedral in dihedrals
        if dihedral.residue.name == "P" && dihedral.dtype == Common.phi
            rotate_dihedral!(state.xyz, dihedral.a2, dihedral.a1, deg2rad(180), Common.phi, dihedral.residue.atoms, dihedral.residue)
            # insert!(dihedral.movable, 1, dihedral.residue.atoms[2])
            # insert!(dihedral.movable, 1, dihedral.residue.atoms[3])
        end
    end
end