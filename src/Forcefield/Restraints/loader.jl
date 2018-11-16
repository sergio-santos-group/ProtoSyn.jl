function load_distance_restraints_from_file(input_file::String, metadata::Common.Metadata; k::Float64 = 1.0)::Vector{DistanceFBR}
    #All distances are in nm.

    #Gather residue information
    resnum2ca::Dict{Int64, Int64} = Dict{Int64, Int64}()
    for residue in metadata.residues
        ca_index = filter(atom -> metadata.atoms[atom].name == "CA", residue.atoms)[1]
        resnum2ca[metadata.atoms[ca_index].res_num] = ca_index
    end

    #Compile restraints
    restraints::Vector{DistanceFBR} = Vector{DistanceFBR}()
    open(input_file) do f
        for line in eachline(f)
            elem = split(line)
            if length(elem) > 0 && elem[1] != "i"
                r1_index = parse(Int64, elem[1])
                r2_index = parse(Int64, elem[2])
                r1 = metadata.residues[r1_index]
                r2 = metadata.residues[r2_index]
                r1_ca = resnum2ca[r1_index]
                r2_ca = resnum2ca[r2_index]
                c = parse(Float64, elem[5])
                if r1.ss == r2.ss && r1.ss == Common.SS.SHEET
                    push!(restraints, DistanceFBR(r1_ca, r2_ca, -Inf, -Inf, 0.4, 0.6, c * k))
                else
                    push!(restraints, DistanceFBR(r1_ca, r2_ca, -Inf, -Inf, 0.5, 1.0, c * k))
                end
            end
        end
    end
    return restraints
end


function compile_dihedral_restraints_from_metadata(metadata::Common.Metadata; k::Float64 = 1.0)::Vector{DihedralFBR}

    eq_angles = Dict(
        Common.SS.SHEET => Dict(Common.DIHEDRAL.phi => deg2rad(-139.0), Common.DIHEDRAL.psi => deg2rad(135.0)),
        Common.SS.HELIX => Dict(Common.DIHEDRAL.phi => deg2rad(-57.0),  Common.DIHEDRAL.psi => deg2rad(-47.0)))

    fr_angle = Dict(Common.SS.SHEET => deg2rad(10), Common.SS.HELIX => deg2rad(20))

    restraints::Vector{DihedralFBR} = Vector{DihedralFBR}()
    for dihedral in metadata.dihedrals
        if dihedral.residue.ss == Common.SS.COIL || dihedral.dtype >= Common.DIHEDRAL.omega
            continue
        end
        r2 = eq_angles[dihedral.residue.ss][dihedral.dtype] - fr_angle[dihedral.residue.ss]
        r3 = eq_angles[dihedral.residue.ss][dihedral.dtype] + fr_angle[dihedral.residue.ss]
        push!(restraints, DihedralFBR(dihedral.a1, dihedral.a2, dihedral.a3, dihedral.a4, -Inf, r2, r3, Inf, k))
    end
    return restraints
end