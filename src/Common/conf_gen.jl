@doc raw"""
    apply_initial_conf!(state::State, dihedrals::Vector{Dihedral})

Apply predefined angles to all dihedrals defined in `dihedrals`, based on the [`Dihedral`](@ref).residue.ss, changing the State.xyz
to apply the secondary structure. The applied angles (in degrees) are the following:

Beta sheet:
PHI = -139.0 | PSI = 135.0
Alpha helix:
PHI = -57.0  | PSI = -47.0

# Examples
```julia-repl
julia> Common.apply_initial_conf(state, dihedrals)
```
"""
function apply_initial_conf!(state::State, dihedrals::Vector{Dihedral})

    #1) Define target values for dihedral angles
    t_angles = Dict(
        beta  => Dict(phi => deg2rad(-139.0), psi => deg2rad(135.0)),
        alpha => Dict(phi => deg2rad(-57.0),  psi => deg2rad(-47.0))
    )

    for dihedral in dihedrals
        dihedral.residue.ss == coil ? continue : nothing
        dihedral.dtype > omega ? continue : nothing

        #2) Calculate current value for the dihedral angle
        current_angle = Aux.calc_dih_angle(
            state.xyz[dihedral.a1, :],
            state.xyz[dihedral.a2, :],
            state.xyz[dihedral.a3, :],
            state.xyz[dihedral.a4, :])

        #3) Calculate necessary displacement to target angles
        displacement = t_angles[dihedral.residue.ss][dihedral.dtype] - current_angle
        
        #4) Apply displacement and rotate
        rotate_dihedral!(state.xyz, dihedral, displacement)
    end
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
            insert!(dihedral.movable, 1, dihedral.residue.atoms[2])
            insert!(dihedral.movable, 1, dihedral.residue.atoms[3])
        end
    end
end