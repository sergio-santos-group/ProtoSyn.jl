@doc raw"""
    load_distance_restraints_from_file(input_file::String, metadata::Common.Metadata[, k::Float64 = 1.0, threshold::Float64 = 0.0])::Vector{DistanceFBR}

Read a contact map file and load the distances and force constants between Cα's in the structure, if the force constant is above the defined `threshold`.
All force constants are multiplied by a factor `k` (Default: 1.0).
By default, the distances between β-sheets are lower than r3=4Å (r4=6Å), while distances between different type of secondary structures or α-helixes are lower than
r3=5Å (r4=10Å).
Return an array of DistanceFBR.

# Examples
```julia-repl
julia> Forcefield.Restraints.load_distance_restraints_from_file(contact_map, metadata, 1e4)
[Forcefield.Restraints.DistanceFBR(a1=1, a2=2, r1=-Inf, r2=-Inf, r3=0.4, r4=0.6, c=1e4),
 Forcefield.Restraints.DistanceFBR(a1=2, a2=3, r1=-Inf, r2=-Inf, r3=0.5, r4=1.0, c=1e4)
 (...)]
```
"""
function load_distance_restraints_from_file(input_file::String, metadata::Common.Metadata; k::Float64 = 1.0, threshold::Float64 = 0.0, min_distance::Float64 = 0.8)::Vector{DistanceFBR}
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
                c = parse(Float64, elem[5])
                if c >= threshold
                    r1_index = parse(Int64, elem[1])
                    r2_index = parse(Int64, elem[2])
                    r1 = metadata.residues[r1_index]
                    r2 = metadata.residues[r2_index]
                    r1_ca = resnum2ca[r1_index]
                    r2_ca = resnum2ca[r2_index]
                    if r1.ss == r2.ss && r1.ss == Common.SS.SHEET
                        push!(restraints, DistanceFBR(r1_ca, r2_ca, -Inf, -Inf, min_distance, min_distance + 0.2, c * k))
                    else
                        push!(restraints, DistanceFBR(r1_ca, r2_ca, -Inf, -Inf, min_distance, min_distance + 0.2, c * k))
                    end
                end
            end
        end
    end
    println("(  PRE) ▲ Loaded $(length(restraints)) contacts")
    return restraints
end


@doc raw"""
    lock_block_bb(metadata::Common.Metadata[, k::Float64 = 1.0, fbw::Float64 = 2.5])::Vector{DihedralFBR}

Add restraints to the PHI and PSI dihedrals of blocks in the structure (as defined in `metadata`).
All force constants are multiplied by a factor `k` (Default: 1.0).
The angles locked are defined in Common.jl and a flat-bottom width (`fbw`) is added on each side of the derised atom. (Default: 10.0 degrees)
Return an array of DihedralFBR.

# Examples
```julia-repl
julia> Forcefield.Restraints.lock_block_bb(metadata, 1e4)
[Forcefield.Restraints.DistanceFBR(a1=1, a2=2, a3=3, a4=4, r1=-180, r2=-175, r3=-165, r4=-160, c=1e4),
 (...)]
```
"""
function lock_block_bb(metadata::Common.Metadata; k::Float64 = 1.0, fbw::Float64 = 2.5)::Vector{DihedralFBR}

    fbw = deg2rad(fbw)/2
    restraints::Vector{DihedralFBR} = Vector{DihedralFBR}()
    for dihd in metadata.dihedrals
        if dihd.residue.ss == Common.SS.COIL || dihd.dtype >= Common.DIHEDRAL.omega
            continue
        end
        r0 = Common.ss2bbd[dihd.residue.ss][dihd.dtype]
        # println(" + DihedralFBR $(dihd.dtype)-$(dihd.residue.name)($(dihd.a1)|$(dihd.a2)|$(dihd.a3)|$(dihd.a4))")
        push!(restraints, DihedralFBR(dihd.a1, dihd.a2, dihd.a3, dihd.a4, r0-(fbw*2), r0-fbw, r0+fbw, r0+(fbw*2), k))
    end
    return restraints
end