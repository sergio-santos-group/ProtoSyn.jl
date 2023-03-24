function measure_chi_similarity(pose::Pose, target_pose::Pose, selection::Opt{AbstractSelection} = nothing, verbose::Bool = true)

    mirror_residues = rn"PHE|TY|ASP"r(pose, gather = true)

    function check_similarity(atom1::Atom, state1::State, atom2::Atom, state2::State)
        atom1_∠ = ProtoSyn.unit_circle(ProtoSyn.getdihedral(state1, atom1))
        atom2_∠ = ProtoSyn.unit_circle(ProtoSyn.getdihedral(state2, atom2))
        d = abs(atom1_∠ - atom2_∠)
        v = d <= 0.6981317007977318 # 40°
        if !v && atom1.container in mirror_residues && atom1.name == ProtoSyn.Peptides.chi_dict[atom1.container.name][end]
            @info "$atom1 = $(rad2deg(atom1_∠)) ≈ $(rad2deg(atom2_∠)) -> $v"
            @info " Attempting mirror ..."
            atom1_∠ -= π
            d = abs(atom1_∠ - atom2_∠)
            v = d <= 0.6981317007977318 # 40°
        end
        @info "$atom1 = $(rad2deg(atom1_∠)) ≈ $(rad2deg(atom2_∠)) -> $v"
        
        return v
    end
    
    pose_residues   = ProtoSyn.promote(selection, Residue)(pose, gather = true)
    target_residues = ProtoSyn.promote(selection, Residue)(pose, gather = true)
    @assert length(pose_residues) === length(target_residues) "The provided selection yielded a different number of selected residues in the pose ($(length(pose_residues))) and the target ($(length(target_residues)))."
    
    same1 = 0
    same2 = 0
    for (pose_residue, target_residue) in zip(pose_residues, target_residues)
        pose_atom_1 = chi"1"(pose_residue, gather = true)
        pose_atom_2 = chi"2"(pose_residue, gather = true)
        target_atom_1 = chi"1"(pose_residue, gather = true)
        target_atom_2 = chi"2"(pose_residue, gather = true)

        if length(pose_atom_1) > 0 && length(target_atom_1) > 0
            global _same1 = check_similarity(pose_atom_1[1], pose.state, target_atom_1[1], target_pose.state)
        end

        if length(pose_atom_2) > 0 && length(target_atom_2) > 0
            global _same2 = check_similarity(pose_atom_2[1], pose.state, target_atom_2[1], target_pose.state)
        end

        if _same1
            same1 += 1
            if _same2
                same2 += 1
            end
        end
    end

    N = length(pose_residues)
    verbose && println("Same Χ1   = $same1/$N")
    verbose && println("Same Χ1+2 = $same2/$N")
    return same1, same2, N
end