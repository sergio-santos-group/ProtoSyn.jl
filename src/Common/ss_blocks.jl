@doc raw"""
    compile_ss_blocks_metadata(metadata::Metadata, ss::String)

Read the provided secondary structure string `ss` and return the [`SecondaryStructureMetadata`](@ref) and [`BlockMetadata`](@ref) from the provided list of residues or dihedrals in `metadata`.

# Examples
```julia-repl
julia> Common.compile_ss_blocks_metadata(metadata, "CCCHHHHCCEEEECCC"
([SecondaryStructureMetadata(ss_type=HELIX, name=HA, I-4 <-> A-7, conf=1),
SecondaryStructureMetadata(ss_type=SHEET, name=BA, A-10 <-> V-13, conf=1)],
[BlockMetadata(atoms=1<->135, pivot=67, range_left=Inf, connector_left=18, connector_right=135),
BlockMetadata(atoms=177<->362, pivot=269, range_left=1.559, connector_left=177, connector_right=362)])
```
"""
function compile_ss_blocks_metadata(metadata::Metadata, ss::String)
    if length(metadata.residues) != 0
        ss     = compile_ss!(metadata.residues, ss)
        blocks = compile_blocks!(metadata.residues, ss)
        return ss, blocks
    elseif length(metadata.dihedrals) != 0
        ss     = compile_ss(metadata.dihedrals, ss)
        blocks = compile_blocks!(metadata.dihedrals, ss)
        return ss, blocks
    else
        error("Tried to compile secondary structure information from both metadata.residues and metadata.dihedrals, but no metadata information was found (Residues = $(length(metadata.residues)) | Dihedrals = $(length(metadata.residues)))")
    end
end

@doc raw"""
    compile_ss_blocks_metadata!(metadata::Metadata, ss::String)

Read the provided secondary structure string `ss` and compile the [`SecondaryStructureMetadata`](@ref) and [`BlockMetadata`](@ref) from the provided list of residues or dihedrals in `metadata`.

# Examples
```julia-repl
julia> Common.compile_ss_blocks_metadata!(metadata, "CCCHHHHCCEEEECCC")
```
See also: [`compile_ss!`](@ref) [`compile_blocks!`](@ref)
"""
function compile_ss_blocks_metadata!(metadata::Metadata, ss::String)
    if length(metadata.residues) != 0
        metadata.blocks = compile_blocks!(metadata.residues, ss)
        return
    elseif length(metadata.dihedrals) != 0
        metadata.blocks = compile_blocks!(metadata.dihedrals, ss)
        return
    else
        error("Tried to compile secondary structure information from both metadata.residues and metadata.dihedrals, but no metadata information was found (Residues = $(length(metadata.residues)) | Dihedrals = $(length(metadata.residues)))")
    end
end

# @doc raw"""
#     compile_ss!(dihedrals::Vector{Dihedral}, ss::String)::Vector{SecondaryStructureMetadata}

# Read the provided secondary structure string `ss` and compile the [`SecondaryStructureMetadata`](@ref) [`Metadata`](@ref) from the provided list of dihedrals.

# # Examples
# ```julia-repl
# julia> Common.compile_ss!(dihedrals, "CCCHHHHCCEEEECCC")
# 2-element Array{ProtoSyn.Common.SecondaryStructureMetadata,1}:
#  SecondaryStructureMetadata(ss_type=HELIX, name=HA, I-4 <-> A-7, conf=1)  
#  SecondaryStructureMetadata(ss_type=SHEET, name=BA, A-10 <-> V-13, conf=1)
# ```
# See also: [`compile_ss_blocks_metadata!`](@ref Common.compile_ss_blocks_metadata!)
# """
# function compile_ss!(dihedrals::Vector{Dihedral}, ss::String)::Vector{SecondaryStructureMetadata}

#     residues = Vector{Residue}()
#     for dihedral in dihedrals
#         if !(dihedral.residue in residues)
#             push!(residues, dihedral.residue)
#         end
#     end
#     return compile_ss!(residues, ss)
# end


# @doc raw"""
#     compile_ss!(residues::Vector{Residue}, ss::String)::Vector{SecondaryStructureMetadata}

# Read the provided secondary structure string `ss` and compile the [`SecondaryStructureMetadata`](@ref) [`Metadata`](@ref) from the provided list of residues.

# # Examples
# ```julia-repl
# julia> Common.compile_ss!(residues, "CCCHHHHCCEEEECCC")
# 2-element Array{ProtoSyn.Common.SecondaryStructureMetadata,1}:
#  SecondaryStructureMetadata(ss_type=HELIX, name=HA, I-4 <-> A-7, conf=1)  
#  SecondaryStructureMetadata(ss_type=SHEET, name=BA, A-10 <-> V-13, conf=1)
# ```
# See also: [`compile_ss_blocks_metadata!`](@ref Common.compile_ss_blocks_metadata!)
# """
# function compile_ss!(residues::Vector{Residue}, ss::String)::Vector{SecondaryStructureMetadata}

#     conv_type = Dict('H' => SS.HELIX, 'E' => SS.SHEET, 'C' => SS.COIL)
#     conv_name = Dict('H' => "HA", 'E' => "BA")

#     sec_str = SecondaryStructureMetadata[]
#     last_ss::Char = ss[1]
#     i_idx::Int64 = 1
#     for (index, curr_ss) in enumerate(ss)
#         if curr_ss != last_ss
#             if last_ss in ['H', 'E']
#                 push!(sec_str, SecondaryStructureMetadata(conv_name[last_ss], residues[i_idx].name, i_idx, residues[index - 1].name, index - 1))
#             end
#             i_idx = index
#         end
#         residues[index].ss = conv_type[curr_ss] # If commented, will not apply ss to residue.ss
#         last_ss = curr_ss
#     end
#     if last_ss in ['H', 'E']
#         push!(sec_str, SecondaryStructureMetadata(conv_name[last_ss], residues[i_idx].name, i_idx, residues[length(ss) - 1].name, length(ss) - 1))
#     end
#     printstyled("(SETUP) ▲ Compiled metadata information of $(length(sec_str)) secondary structures\n", color = 9)
#     return sec_str
# end

@doc raw"""
    compile_blocks!(dihedrals::Vector{Dihedral}, ss::String)::Vector{BlockMetadata}

Read the provided secondary structure string `ss` and compile the [`BlockMetadata`](@ref) [`Metadata`](@ref) from the provided list of dihedrals.

# Examples
```julia-repl
julia> Common.compile_blocks!(dihedrals, "HHHHCCEEEECCC")
2-element Array{ProtoSyn.Common.BlockMetadata,1}:
 BlockMetadata(atoms=1<->135, pivot=67, range_left=Inf, connector_left=18, connector_right=135)
 BlockMetadata(atoms=177<->362, pivot=269, range_left=1.559, connector_left=177, connector_right=362)
```
See also: [`compile_ss_blocks_metadata!`](@ref Common.compile_ss_blocks_metadata!)
"""
function compile_blocks!(dihedrals::Vector{Dihedral}, ss::String)::Vector{BlockMetadata}

    residues = Vector{Residue}()
    for dihedral in dihedrals
        if !(dihedral.residue in residues)
            push!(residues, dihedral.residue)
        end
    end
    return compile_blocks!(residues, ss)
end

# Takes information from string and puts the correct SS.TYPE on each residue.ss
# TODO: Documentation
function apply_ss_to_residues!(residues::Vector{Residue}, ss::String)

    conv = Dict('H' => SS.HELIX, 'E' => SS.SHEET, 'C' => SS.COIL)

    for (index, s) in enumerate(ss)
        residues[index].ss = conv[s]
    end
end

# Takes information from string, puts the correct SS.TYPE on each
# residue.ss and returns the NON-COIL blocks
# TODO: Documentation
function compile_blocks!(residues::Vector{Residue}, ss::String)::Vector{BlockMetadata}
    d1 = 0.24510
    d2 = 0.11536

    apply_ss_to_residues!(residues, ss)
    ss_array = Common.retrieve_ss_from_residues(residues)

    blocks = Vector{BlockMetadata}()
    for (ss_index, ss) in enumerate(ss_array)
        if ss.type != SS.COIL

            atoms = Vector{Int64}()
            for residue in residues[ss.i_res_num:ss.f_res_num]
                atoms = vcat(atoms, residue.atoms)
            end

            # Connector_left / Connector_right = atom ID used to
            # calculate the current distance of the left coil.
            # This distance is calculated as the position of
            # Connector_left - Connector_right (of the previous block to the left)
            # This distance is compared to `range_left`: Range_left is the
            # maximum allowed distance that the left coil can support when
            # fully extended. 
            # Range left must be in nm.
            c_left  = atoms[1]
            c_right = atoms[end - 1]
            coil_left = ss_array[ss_index - 1]
            dn = coil_left.f_res_num - coil_left.i_res_num + 1
            range_left = (d1 * dn + d2 * (dn +1))
                        
            # Special cases (First and last block take the residual coils on
            # each side as part of the block for rotations, etc)
            if ss_index == 2 && ss_array[1].type == SS.COIL
                residual_atoms = Vector{Int64}()
                for residue in residues[1:(ss.i_res_num - 1)]
                    residual_atoms = vcat(residual_atoms, residue.atoms)
                end
                atoms = vcat(residual_atoms, atoms)
            end
            if ss_index == (length(ss_array) - 1) && ss_array[end].type == SS.COIL
                for residue in residues[(ss.f_res_num + 1):end]
                    atoms = vcat(atoms, residue.atoms)
                end
            end

            push!(blocks, BlockMetadata(ss.type, atoms, range_left, c_left, c_right))
        end
    end
    printstyled("(SETUP) ▲ Compiled metadata information of $(length(blocks)) blocks\n", color = 9)
    return blocks
end

# Scans the molecule residues and returns and the secondary structures including
# coils. Used for Printing.
function retrieve_ss_from_residues(residues::Vector{Residue})::Vector{SecondaryStructureMetadata}
    
    ss_array = Vector{SecondaryStructureMetadata}()
    
    index = 1
    c = residues[index]

    while true
        ss = SecondaryStructureMetadata(c.ss, c.name, index, c.name, index)
        while true
            index += 1
            (c.next == nothing || c.next.ss != c.ss) && break
            ss.f_res_name = c.next.name
            ss.f_res_num  = index
            c = c.next
        end
        push!(ss_array, ss)
        c = c.next
        c == nothing && break
    end
    return ss_array
end