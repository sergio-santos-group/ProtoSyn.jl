# ----------------------------------------------------------------------------------------------------------
#                                                 LOADERS

# #TODO: Verify function
# @doc raw"""
#     load_from_gro(i_file::String)::Common.State

# Return a new [`Common.State`](@ref) by loading the atom positions and metadata from the input .gro file.
# As a default, `state.energy` is [`NullEnergy`](@ref) and `state.forces` are set to zero.
# If `compile_metadata` flag is set to true, the returned Metadata object is compiled from the existing  information in the GRO file.

# # Examples
# ```julia-repl
# julia> Common.load_from_gro("molecule.gro")
# Common.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], metadata=(...))
# ```
# See also: [`load_from_pdb`](@ref)
# """
# function load_from_gro(i_file::String, compile_metadata = true)::Tuple{State, Metadata}

#     #Initialize empty arrays
#     xyz   = Vector{Array{Float64, 2}}()
#     if compile_metadata atoms = Vector{AtomMetadata}() end

#     #Read file (from file name) and recover XYZ and ATOM NAMES
#     open(i_file, "r") do f
#         for (index, line) in enumerate(eachline(f))
#             elem = split(line)
#             if length(elem) > 3 && index > 2
#                 push!(xyz, map(x -> parse(Float64, x), [elem[4] elem[5] elem[6]]))
#                 if compile_metadata
#                     push!(atoms, AtomMetadata(index=index, name=string(elem[2]), elem=string(elem[2]), res_num=parse(Int64, strip(line[1:5])), res_name=string(strip(line[6:8]))))
#                 end
#             end
#         end
#     end

#     n = length(xyz)
#     if compile_metadata
#         metadata = Metadata(atoms = atoms)
#     else
#         metadata = Metadata()
#     end
#     return State(n, Common.Energy(), vcat(xyz...), zeros(n, 3)), metadata
# end


@doc raw"""
    load_from_pdb(i_file::String[, compile_metadata::Bool = true])::Common.State

Return a new [`Common.State`](@ref) by loading the atom positions and metadata from the input .pdb file.
As a default, `state.energy` is [`NullEnergy`](@ref) and `state.forces` are set to zero.
If `compile_metadata` flag is set to true, the returned Metadata object is compiled from the existing  information in the PDB file.

# Examples
```julia-repl
julia> Common.load_from_pdb("molecule.pdb")
Common.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], metadata=(...))
```
See also: [`load_from_gro`](@ref)
"""
function load_from_pdb(i_file::String, compile_metadata = true)::Tuple{State, Metadata}

    xyz = Vector{Vector{Float64}}()
    md = Metadata()
    open(i_file, "r") do f
        for line in eachline(f)
            if startswith(line, "ATOM")
                push!(xyz, [0.1*parse(Float64, line[n:n+7]) for n=31:8:47])
            end
        end
        if compile_metadata
            md = load_metadata(f)
        end
    end
    
    state = State(size(xyz, 1))
    state.xyz = reshape(vcat(xyz...),:,3)
    # state.xyz = collect(reshape(vcat(xyz...),:,3)')

    state, md
    # if compile_metadata
    #     metadata = Metadata(atoms = atoms)
    #     compile_residue_metadata!(metadata)
    #     compile_dihedral_metadata!(metadata)
    # else
    #     metadata = Metadata()
    # end

    # return State(n, Energy(), vcat(xyz...), zeros(n, 3), NonBondedList(n)), metadata
    #n = length(xyz)
    # return State(n, Energy(), vcat(xyz...), zeros(n, 3), nothing), metadata
    #return State(n, Energy(), reshape(vcat(xyz...),3,:), zeros(n, 3), nothing), metadata
end

function load_metadata(io::IOStream)
    cur_res_id = -1
    atoms = Vector{AtomMetadata}()
    residues = Vector{Residue}()

    seekstart(io)
    for line in eachline(io)
        if startswith(line, "ATOM")
            rindex = parse(Int, line[23:26])
            rname = string(strip(line[18:20]))
            if cur_res_id != rindex
                cur_res_id = rindex
                residue = Residue(index=rindex, name=rname)
                push!(residues, residue)
            end
            amd = AtomMetadata(
                index  = parse(Int, strip(line[6:11])),
                name   = string(strip(line[14:16])),
                symbol = string(strip(line[77:78])),
                residue = residues[end]
            )
            push!(atoms, amd)
        elseif startswith(line, "CONECT")
            indices = map(s->parse(Int, s), [line[n:n+4] for n=7:5:length(line)])
            atoms[indices[1]].connects = indices[2:end]
        end
    end

    Metadata(atoms=atoms, residues=residues)
end

# @doc raw"""
#     load_metadata_from_json(json_file::String)

# Parse a JSON file containing the dihedral and residue topology. Return a [`Dihedral`](@ref) array
# and a [`Residue`](@ref) array.

# # Examples
# ```julia-repl
# julia> Mutators.Diehdral.load_topology("mc.json")
# (ProtoSyn.Mutators.Dihedral.NewDihedral[...], ProtoSyn.Common.Residue[...])
# ```
# See also: [`Aux.read_JSON`](@ref)
# """
# function load_metadata_from_json(json_file::String)

#     p::Dict{String, Any} = Aux.read_JSON(json_file)
#     tmp_residues = Dict(d["n"] => Residue(convert(Vector{Int64}, d["atoms"]), d["next"], d["type"]) for d in p["residues"])
    
#     str2enum = Dict(string(s) => s for s in instances(DIHEDRAL.TYPE))
    
#     dihedrals = [
#         Dihedral(d["a1"], d["a2"], d["a3"], d["a4"],
#             d["movable"], tmp_residues[d["parent"]], str2enum[lowercase(d["type"])])
#         for d in p["dihedrals"]
#     ]
    
#     # Set correct references for dihedrals next
#     for residue in values(tmp_residues)
#         residue.next = get(tmp_residues, residue.next, nothing)
#     end


#     #Export residues as an ordered dictionary
#     residues = Vector{Residue}()
#     for i in 1:length(tmp_residues)
#         push!(residues, tmp_residues[i])
#     end

#     return dihedrals, residues
# end


# @doc raw"""
#     compile_residue_metadata!(metadata::Metadata)

# Add residue information to `metadata`, compiling the avaliable information from the metadata (requires pre-existing information regarding the AtomMetadata).

# # Examples
# ```julia-repl
# julia> Common.compile_residue_metadata!(metadata)
# ```
# See also: [`load_metadata_from_json`](@ref)
# """
# function compile_residue_metadata!(metadata::Metadata)

#     # Verify input
#     if length(metadata.atoms) <= 0
#         error("Metadata needs to have AtomMetadata information in order to compile residues.")
#     end

#     # This function will iterate over each atom metadata in reverse order. This allows to define the `residue.next` as the previously defined residue in the iteration,
#     # since it is being looked up in reverse. The first residue defined does not have a `next`, so it is exeptionally set to `nothing`. A new residue is identified
#     # whenever the atom being iterated over has a different `res_num` than the previous one. Newly identified residues are appended to the begining of the list, so
#     # they are returned in order.
#     residue_metadata = Vector{Residue}()
    
#     # Iterate over the atom list in reverse order
#     residues = Dict{Int64, Vector{AtomMetadata}}()
#     for atom in metadata.atoms
#         if !(atom.res_num in keys(residues))
#             residues[atom.res_num] = Vector{AtomMetadata}()
#         end
#         push!(residues[atom.res_num], atom)
#     end

#     # Iterate over the reversed residues list
#     for key in reverse(sort(collect(keys(residues))))
#         residue = residues[key]
#         sort!(residue, by = atom -> atom.index)
#         if length(residue_metadata) == 0
#             # Exceptionally, the "first" residue (actually last, since it's in reverse order) will have `nothing` as its `next`.
#             next = nothing
#         else
#             # All other residues have the previously defined residue as its `next`.
#             next = residue_metadata[1]
#         end
#         # Add this residue to the residues list, appending it to the begining of the list.
#         res = Residue([atom.index for atom in residue], next, residue[1].res_name, SS.COIL)
#         insert!(residue_metadata, 1, res)
#         for atom in residue
#             atom.residue = residue_metadata[1]
#         end
#     end
#     metadata.residues = residue_metadata
# end


# @doc raw"""
#     compile_dihedral_metadata!(metadata::Metadata)

# Add dihedrals information to `metadata`, compiling the avaliable information from the metadata (requires pre-existing information regarding the AtomMetadata).

# # Examples
# ```julia-repl
# julia> Common.compile_dihedral_metadata!(metadata)
# ```
# See also: [`load_metadata_from_json`](@ref)
# """
# function compile_dihedral_metadata!(metadata::Metadata)

#     function find_intra_residue_movables!(atoms::Vector{AtomMetadata}, a3::Int64, current_res::Int64, current_list::Vector{Int64}, exclude::Int64)::Vector{Int64}
#         # This function receives the full atoms metadata list so it can look up the current atom `res_num` and make sure the recursive search does not extend to other
#         # residues AND to search for the `connects` records. The recursive search is also stopped if the connect is already in the `movables` list or is equal to the
#         # exclude atom (that identifies the previously looked up atom in the chain).
#         push!(current_list, a3)
#         if atoms[a3].connects != nothing
#             for connect in atoms[a3].connects
#                 if (connect in current_list || atoms[connect].res_num != current_res || connect == exclude)
#                     continue
#                 end
#                 find_intra_residue_movables!(atoms, connect, current_res, current_list, a3)
#             end
#         end
#         return current_list
#     end

#     #Verify input
#     if length(metadata.atoms) <= 0
#         error("Metadata needs to have AtomMetadata information in order to compile dihedrals.")
#     end

#     # This function iterates over each residue (ignoring Proline) and identifies the side-chain dihedrals. The atoms list (containing each atom metadata)
#     # is necessary to identify atom's names.
#     dihedrals = Vector{Dihedral}()
    
#     # Identify residues from atom data
#     residues = Dict{Int64, Vector{AtomMetadata}}()
#     for atom in metadata.atoms
#         if !(atom.res_num in keys(residues))
#             residues[atom.res_num] = Vector{AtomMetadata}()
#         end
#         push!(residues[atom.res_num], atom)
#     end

#     for residue in sort(collect(values(residues)), by = atoms -> atoms[1].index)
        
#         # Create a new atom_name -> atom conversion dictionary for each residue iterated over.
#         name2atom = Dict(atom.name => atom for atom in residue)

#         # Identify Backbone Dihedrals
#         n  = name2atom["N"]
#         ca = name2atom["CA"]
#         c  = name2atom["C"]
#         # PHI
#         # if n.connects != nothing && residue[1].res_name != "PRO"
#         if n.connects != nothing
#             prev_c = filter(connect -> metadata.atoms[connect].name == "C", n.connects)
#             if length(prev_c) > 0
#                 movables = find_intra_residue_movables!(metadata.atoms, ca.index, ca.res_num, Vector{Int64}(), n.index)
#                 push!(dihedrals, Dihedral(prev_c[1], n.index, ca.index, c.index, sort(movables), ca.residue, DIHEDRAL.phi))
#             end
#         end
#         # PSI
#         if c.connects != nothing
#             next_n = filter(connect -> metadata.atoms[connect].name == "N", c.connects)
#             if length(next_n) > 0
#                 movables = find_intra_residue_movables!(metadata.atoms, c.index, c.res_num, Vector{Int64}(), ca.index)
#                 push!(dihedrals, Dihedral(n.index, ca.index, c.index, next_n[1], sort(movables), c.residue, DIHEDRAL.psi))
#                 # OMEGA
#                 next_n = metadata.atoms[next_n[1]]
#                 if next_n.connects != nothing
#                     next_ca = filter(connect -> metadata.atoms[connect].name == "CA", next_n.connects)
#                     if length(next_ca) > 0
#                         # Omega dihedrals don't have movables
#                         push!(dihedrals, Dihedral(ca.index, c.index, next_n.index, next_ca[1], Vector{Int64}(), c.residue, DIHEDRAL.omega))
#                     end
#                 end
#             end
#         end

#         # Identify Side-chain Dihedrals
#         # Ignore Proline
#         if residue[1].res_name == "PRO"
#             continue
#         end
#         # Start the dihedral path
#         path = ["N", "CA"]
#         # Based on the `residue.name`, the possible path taken by the iteration is different. This accounts for cyclic aminoacids.
#         possible_queries = residue[1].res_name in ("HIS", "PHE", "TRP", "TYR") ? "BGD" : "BGDEZH"
#         # Join all atom names in this residue in a single string, dividing atom names with ":". This will allow identification of the dihedral using regular expressions.
#         atnames = string(":", join(map(atom -> atom.name, residue), ":"), ":")
#         # Iterate over the given possible path for the side-chain.
#         for (chi_number, query) in enumerate(possible_queries)
#             # The regular expression tries to identify the current query in the list of atom names. This is an overpowered version of looking for a key in the dictionary
#             # since it account for possible numbers in the atom name and different elements (C vs N for example).
#             m = match(Regex("(.*):(?P<name>[A-GI-Z]$query[1]?):(.*)"), atnames)
#             # If the current query is found, add this atom name to the dihedral path being built.
#             if m != nothing
#                 push!(path, (m[:name]))
#             end
#             # If the dihedral path reaches 4 atoms (meaning a dihedral was identified), add a new Dihedral to the dihedrals list.
#             if length(path) == 4
#                 # The previously defined atom_name -> atom_index conversion dictionary is employed to define the dihedral by its *global* indices.
#                 a1, a2, a3, a4 = map(atom_name -> name2atom[atom_name].index, path)
#                 # The dihedral movables are identified recursively using find_intra_residue_movables!() function.
#                 movables = find_intra_residue_movables!(metadata.atoms, a3, metadata.atoms[a3].res_num, Vector{Int64}(), a2)
#                 # The dihedral type is defined in order, starting in chi1, chi2, etc. Since the `enumerate` function, in Julia, starts in 1 and the first dihedral will
#                 # only be identified in the second iteration of the loop (since on the first iteration only 3 atoms of the path have been identified), `chi_number - 1`
#                 # is used.
#                 push!(dihedrals, Dihedral(a1, a2, a3, a4, sort(movables), residue[1].residue, DIHEDRAL.TYPE(chi_number - 1)))
#                 # Delete the begining of the dihedral path. If a new atom is identified in the side-chain, the a2 of this dihedral becomes a1, etc
#                 deleteat!(path, 1)
#             end
#         end
#     end
#     metadata.dihedrals = dihedrals
# end

# @doc raw"""
#     compile_sidechains_metadata(metadata::Metadata, rl::Dict{String, Any})

# Return sidechains information, compiling the avaliable information from the [`Metadata`](@ref)
# (requires pre-existing information regarding the Rotamer Library, from an external JSON file, where rotamer angles are in degrees).

# # Examples
# ```julia-repl
# julia> Common.compile_sidechains_metadata!(metadata, aux.read_JSON(rotamer_library))
# 45-element Array{ProtoSyn.Common.SidechainMetadata,1}:
#  [SidechainMetadata(dihedrals=3, rotamers=2, weights=[0.0254237, 0.0360169],
#  (...)]
# ```
# """
# function compile_sidechains_metadata(metadata::Metadata, rl::Dict{String, Any})::Vector{SidechainMetadata}
#     if length(metadata.dihedrals) != 0
#         return compile_sidechains_metadata(metadata.dihedrals, rl)
#     else
#         error("Tried to compile sidechain information from metadata.dihedrals (Length: $(length(metadata.dihedrals))), but no dihedral information was found.")
#     end
# end

# @doc raw"""
#     compile_sidechains_metadata!(metadata::Metadata, rl::Dict{String, Any})

# Add sidechains information to `metadata`, compiling the avaliable information from the metadata
# (requires pre-existing information regarding the Rotamer Library, from an external JSON file, where rotamer angles are in degrees).

# # Examples
# ```julia-repl
# julia> Common.compile_sidechains_metadata!(metadata, aux.read_JSON(rotamer_library))
# ```
# """
# function compile_sidechains_metadata!(metadata::Metadata, rl::Dict{String, Any})
#     if length(metadata.dihedrals) != 0
#         metadata.sidechains  = compile_sidechains_metadata(metadata.dihedrals, rl)
#     else
#         error("Tried to compile sidechain information from metadata.dihedrals (Length: $(length(metadata.dihedrals))), but no dihedral information was found.")
#     end
# end

# @doc raw"""
#     compile_sidechains_metadata(dihedrals::Vecotr{Dihedral}, rl::Dict{String, Any})

# Return sidecahin information, as a vector of [`Dihedral`](@ref) vectors, compiling the avaliable information from the provided [`Dihedral`](@ref)s list

# # Examples
# ```julia-repl
# julia> Common.compile_sidechains_metadata!(metadata.dihedrals)
# 45-element Array{Array{ProtoSyn.Common.Dihedral,1}, 1}:
#  [Dihedral(a1=2, a2=3, a3=5, a4=7, movable=[5, 6], residue=Common.Residue(atoms=[1, 2, 3, 4, 5, 6], next=V, name=A), type=phi),
#  (...)]
# ```
# See also: [`SidechainMutator`](@ref Sidechain)
# """
# function compile_sidechains(dihedrals::Vector{Dihedral})::Vector{Vector{Dihedral}}
#     if length(dihedrals) == 0
#         error("Tried to compile sidechain information from the provided dihedrals list (Length: $(length(dihedrals))), but no dihedral information was found.")
#     end

#     all_sidechain_dihedrals = filter(x -> x.dtype > Common.DIHEDRAL.omega, dihedrals)
#     cur_res                 = all_sidechain_dihedrals[1].residue
#     cur_res_index           = 1
#     sidechain_dihedrals     = Vector{Dihedral}()
#     sidechains              = Vector{Vector{Dihedral}}()

#     for dihedral in all_sidechain_dihedrals
#         if cur_res == dihedral.residue
#             push!(sidechain_dihedrals, dihedral)
#         else
#             push!(sidechains, sidechain_dihedrals)
#             sidechain_dihedrals = Vector{Dihedral}()
#             push!(sidechain_dihedrals, dihedral)
#             cur_res = dihedral.residue
#         end
#     end
#     printstyled("(SETUP) ▲ Compiled metadata information of $(length(sidechains)) sidechains\n", color = 9)
#     return sidechains
# end