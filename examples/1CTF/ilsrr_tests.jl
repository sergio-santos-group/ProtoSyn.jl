function calculate_dihedral(a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64}, a4::Vector{Float64})::Float64

    v12 = a2 - a1
    v23 = a3 - a2
    v34 = a4 - a3
    v123 = cross(v12, v23)
    v234 = cross(v23, v34)
    return atan(dot(cross(v123, v234), v23)/sqrt(dot(v23, v23)), dot(v123, v234))
end

function get_FBR_from_dihedral(dihedral::Common.Dihedral, restraint_list::Vector{Forcefield.Restraints.DihedralFBR})::Union{Nothing, Forcefield.Restraints.DihedralFBR}

    for restraint in restraint_list
        if (restraint.a1 == dihedral.a1) && (restraint.a2 == dihedral.a2) && (restraint.a3 == dihedral.a3) && (restraint.a4 == dihedral.a4)
            return restraint
        end
    end
    return nothing
end

function identify_angles_under_stress(dihedrals::Vector{Common.Dihedral}, restraints::Vector{Forcefield.Restraints.DihedralFBR}, threshold::Float64 = 1.0)

    count::Int64 = 0
    for dihedral in dihedrals
        fbr = get_FBR_from_dihedral(dihedral, dihedral_restraints)
        if dihedral.dtype <= Common.DIHEDRAL.omega && fbr != nothing
            angle = calculate_dihedral(state.xyz[dihedral.a1, :], state.xyz[dihedral.a2, :], state.xyz[dihedral.a3, :], state.xyz[dihedral.a4, :])
            if angle < fbr.r2 - deg2rad(threshold) || angle > fbr.r3 + deg2rad(threshold)
                count += 1
                println(@sprintf "%s(%3d-%3d-%3d-%3d)%s -> %7.2f (%7.2f|%7.2f|%7.2f|%7.2f)" dihedral.dtype dihedral.a1 dihedral.a2 dihedral.a3 dihedral.a4 dihedral.residue.name rad2deg(angle) rad2deg(fbr.r1) rad2deg(fbr.r2) rad2deg(fbr.r3) rad2deg(fbr.r4))
            end
        end
    end
    if count == 0
        println("No dihedrals were identified as being under restrictional stress.")
    else
        println("$count dihedrals were identified as being under restrictional stress.")
    end
end