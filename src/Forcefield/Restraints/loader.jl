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


function lock_block_bb(metadata::Common.Metadata; k::Float64 = 1.0, fbw::Float64 = 10.0)::Vector{DihedralFBR}

    fbw = deg2rad(fbw)/2
    restraints::Vector{DihedralFBR} = Vector{DihedralFBR}()
    for dihd in metadata.dihedrals
        if dihd.residue.ss == Common.SS.COIL || dihd.dtype >= Common.DIHEDRAL.omega
            continue
        end
        r0 = Common.ss2bbd[dihd.residue.ss][dihd.dtype]
        push!(restraints, DihedralFBR(dihd.a1, dihd.a2, dihd.a3, dihd.a4, r0-(fbw*2), r0-fbw, r0+fbw, r0+(fbw*2), k))
    end
    return restraints
end