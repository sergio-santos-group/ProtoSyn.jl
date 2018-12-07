@doc raw"""
    compile_ss(dihedrals::Vector{Dihedral}, ss::String)::Vector{SecondaryStructureMetadata}
"""
function compile_ss(dihedrals::Vector{Dihedral}, ss::String)::Vector{SecondaryStructureMetadata}

    residues = Vector{Residue}()
    for dihedral in dihedrals
        if !(dihedral.residue in residues)
            push!(residues, dihedral.residue)
        end
    end
    return compile_ss(residues, ss)
end


@doc raw"""
    compile_ss(residues::Vector{Residue}, ss::String)::Vector{SecondaryStructureMetadata}

Read the provided secondary structure string `ss` and compile the [`SecondaryStructureMetadata`](@ref) [`Metadata`](@ref) from the provided list of residues/dihedrals.

# Examples
```julia-repl
julia> Common.compile_ss(dihedrals, "CCCHHHHCCEEEECCC")
2-element Array{ProtoSyn.Common.SecondaryStructureMetadata,1}:
 SecondaryStructureMetadata(ss_type=HELIX, name=HA, I-4 <-> A-7, conf=1)  
 SecondaryStructureMetadata(ss_type=SHEET, name=BA, A-10 <-> V-13, conf=1)
```
"""
function compile_ss(residues::Vector{Residue}, ss::String)::Vector{SecondaryStructureMetadata}

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
        residues[index].ss = conv_type[curr_ss] # If commented, will not apply ss to residue.ss
        last_ss = curr_ss
    end
    println("(  PRE) ▲ Compiled metadata information of $(length(sec_str)) secondary structures")
    return sec_str
end


function compile_blocks(residues::Vector{Residue}, ss::String)::Vector{BlockMetadata}

    function count_blocks(string::String)::Int64
        count::Int64 = 0
        curr::Char = string[1]
        for l in string
            if l != curr && l in ['H', 'E']
                count += 1
            end
            curr = l
        end
        return count
    end

    n_blocks = count_blocks(ss)
    blocks = BlockMetadata[]
    last_ss::Char = ss[1]
    atoms = Int64[]
    connector_left::Int64 = -1
    i_idx::Int64 = 1
    range_left::Float64 = -1.0
    # single_residue_range::Float64 = 3.636
    connection_range::Float64 = 1.2
    for (index, curr_ss) in enumerate(ss)
        if curr_ss != last_ss
            if last_ss in ['H', 'E']
                connector_right = residues[index - 1].atoms[length(residues[index - 1].atoms)]
                pivot = atoms[floor(Int64, length(atoms)/2)]
                if length(blocks) == n_blocks - 1
                    atoms = vcat(atoms, residues[index].atoms)
                end
                push!(blocks, BlockMetadata(atoms, pivot, range_left, connector_left, connector_right))
                i_idx = index
            else
                connector_left = residues[index].atoms[1]
                if length(blocks) > 0
                    range_left = floor((((index - i_idx) * (connection_range * 3)) + connection_range) / 10, digits = 3)
                else
                    range_left = Inf
                end
            end
            if length(blocks) > 0
                atoms = Int64[]
            end
        end
        last_ss = curr_ss
        atoms = vcat(atoms, residues[index].atoms)
    end
    println("(  PRE) ▲ Compiled metadata information of $(length(blocks)) blocks")
    return blocks
end