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