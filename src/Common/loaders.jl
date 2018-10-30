# ----------------------------------------------------------------------------------------------------------
#                                                 LOADERS

#TODO: Verify function
@doc raw"""
    load_from_gro(i_file::String)::Common.State

    Return a new [`Common.State`](@ref) by loading the atom positions and metadata from the input .gro file.
As a default, `state.energy` is [`NullEnergy`](@ref) and `state.forces` are set to zero.

# Examples
```julia-repl
julia> Common.load_from_gro("molecule.gro")
Common.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], metadata=(...))
```
See also: [`load_from_pdb`](@ref)
"""
function load_from_gro(i_file::String)::Common.State

    #Initialize empty arrays
    xyz     = Vector{Array{Float64, 2}}()
    metadata = Vector{AtomMetadata}()

    #Read file (from file name) and recover XYZ and ATOM NAMES
    open(i_file, "r") do f
        for (index, line) in enumerate(eachline(f))
            elem = split(line)
            if length(elem) > 3 && index > 2
                push!(xyz, map(x -> parse(Float64, x), [elem[4] elem[5] elem[6]]))
                push!(metadata, AtomMetadata(string(elem[2]), elem=string(elem[2]), res_num=parse(Int64, strip(line[1:5])), res_name=string(strip(line[6:8]))))
            end
        end
    end

    n = length(xyz)
    return Common.State(n, Common.Energy(), vcat(xyz...), zeros(n, 3), Metadata(metadata))
end


@doc raw"""
    load_from_pdb(i_file::String)::Common.State

Return a new [`Common.State`](@ref) by loading the atom positions and metadata from the input .pdb file.
As a default, `state.energy` is [`NullEnergy`](@ref) and `state.forces` are set to zero.

# Examples
```julia-repl
julia> Common.load_from_pdb("molecule.pdb")
Common.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], metadata=(...))
```
See also: [`load_from_gro`](@ref)
"""
function load_from_pdb(i_file::String)

    xyz      = Vector{Array{Float64, 2}}()
    metadata = Vector{AtomMetadata}()

    open(i_file, "r") do f
        for line in eachline(f)
            # if length(line) > 6 && line[1:6] == "ATOM  "
            if startswith(line, "ATOM")
                push!(xyz, map(x -> 0.1*parse(Float64, x), [line[31:38] line[39:46] line[47:54]]))
                push!(metadata, AtomMetadata(string(strip(line[14:16])),
                    elem = string(strip(line[77:78])),
                    res_num = parse(Int64, line[23:26]),
                    res_name = string(strip(line[18:20])),
                    chain_id = string(line[22]),
                    connects = nothing))
            elseif startswith(line, "CONECT")
                elem = split(line)
                metadata[parse(Int64, elem[2])].connects = map(x -> parse(Int64, x), elem[3:end])
            end
        end
    end

    n = length(xyz)
    return State(n, Energy(), vcat(xyz...), zeros(n, 3)), Metadata(metadata)
end


@doc raw"""
    load_topology(p::Dict{String, Any})

Parse a dictionary containing the dihedral and residue topology. Return a [`Dihedral`](@ref) array
and a [`Residue`](@ref) array.

# Examples
```julia-repl
julia> Mutators.Diehdral.load_topology(p)
(ProtoSyn.Mutators.Dihedral.NewDihedral[...], ProtoSyn.Common.Residue[...])
```
See also: [`Aux.read_JSON`](@ref)
"""
function load_topology(p::Dict{String, Any})

    tmp_residues = Dict(d["n"] => Residue(convert(Vector{Int64}, d["atoms"]), d["next"], d["type"]) for d in p["residues"])
    
    str2enum = Dict(string(s) => s for s in instances(DIHEDRALTYPE))
    
    dihedrals = [
        Dihedral(d["a1"], d["a2"], d["a3"], d["a4"],
            d["movable"], tmp_residues[d["parent"]], str2enum[lowercase(d["type"])])
        for d in p["dihedrals"]
    ]
    
    # Set correct references for dihedrals next
    for residue in values(tmp_residues)
        residue.next = get(tmp_residues, residue.next, nothing)
    end


    #Export residues as an ordered dictionary
    residues = Vector{Residue}()
    for i in 1:length(tmp_residues)
        push!(residues, tmp_residues[i])
    end

    return dihedrals, residues
end


@doc raw"""
    compile_topology_from_metadata(atoms::Vector{AtomMetadata})::Tuple{Vector{Residue}, Vector{Dihedral}}
    compile_topology_from_metadata(metadata::Metadata)::Tuple{Vector{Residue}, Vector{Dihedral}}

Return both a residue and dihedrals list, compiling the avaliable information from the metadata.    

# Examples
```julia-repl
julia> residues, dihedrals = Common.compile_topology_from_metadata(metadata.atoms)
(...)

julia> residues, dihedrals = Common.compile_topology_from_metadata(metadata)
(...)
```
See also: [`load_topology`](@ref)
"""
function compile_topology_from_metadata(atoms::Vector{AtomMetadata})::Tuple{Vector{Residue}, Vector{Dihedral}}

    # This function returns the residues and dihedrals list of a molecule by parsing the atom motadata identified from the given PDB file.
    # As a benchmark, identification of both residues and dihedrals in the 1CTF metadata is done in ~0.33 seconds.
    residues = compile_residues_from_metadata(atoms)
    bb_dihedrals = compile_backbone_dihedrals_from_metadata(atoms, residues)
    sc_dihedrals = compile_sidechain_dihedrals_from_metadata(atoms, residues)
    
    #Since the mechanisms to identify backbone dihedrals and side-chain dihedrals are different, the two come separated in the final dihedrals list.
    return residues, vcat(bb_dihedrals, sc_dihedrals)
end
compile_topology_from_metadata(metadata::Metadata)::Tuple{Vector{Residue}, Vector{Dihedral}} = compile_topology_from_metadata(metadata.atoms)


@doc raw"""
    compile_residues_from_metadata(atoms::Vector{AtomMetadata})::Vector{Residue}
    compile_residues_from_metadata(metadata::Metadata)::Vector{Residue}

Return a residue list, compiling the avaliable information from the metadata.

# Examples
```julia-repl
julia> residues = Common.compile_residues_from_metadata(metadata.atoms)
(...)

julia> residues = Common.compile_residues_from_metadata(metadata)
(...)
```
See also: [`load_topology`](@ref)
"""
function compile_residues_from_metadata(atoms::Vector{AtomMetadata})::Vector{Residue}

    # This function will iterate over each atom metadata in reverse order. This allows to define the `residue.next` as the previously defined residue in the iteration,
    # since it is being looked up in reverse. The first residue defined does not have a `next`, so it is exeptionally set to `nothing`. A new residue is identified
    # whenever the atom being iterated over has a different `res_num` than the previous one. Newly identified residues are appended to the begining of the list, so
    # they are returned in order.
    residues = Vector{Residue}()
    residue_atoms = Vector{Int64}()
    
    # Iterate over the atom list in reverse order
    curr_res_num = atoms[length(atoms)].res_num
    for index in length(atoms):-1:1

        # A new residue is identified when the `atom.res_num` is different than the previous one.
        if atoms[index].res_num != curr_res_num
            # Exceptionally, the "first" residue (actually last, since it's in reverse order) will have `nothing` as its `next`.
            if length(residues) == 0
                insert!(residues, 1, Residue(reverse(residue_atoms), nothing, atoms[index + 1].res_name, SS.COIL))
            else
                # Add this residue to the residues list, appending it to the begining of the list.
                # All other residues have the previously defined residue as its `next`.
                insert!(residues, 1, Residue(reverse(residue_atoms), residues[1], atoms[index + 1].res_name, SS.COIL))
            end
            # Redefine the current residue number being looked up and create a new empty list of residue atoms.
            curr_res_num = atoms[index].res_num
            residue_atoms = Vector{Int64}()
        end
        # Add the current atom to the residue atoms' list.
        push!(residue_atoms, index)
    end
    # Append the first residue to the begining of the list
    insert!(residues, 1, Residue(reverse(residue_atoms), residues[1], atoms[1].res_name, SS.COIL))
    return residues
end
compile_residues_from_metadata(metadata::Metadata)::Vector{Residue} = compile_residues_from_metadata(metadata.atoms)


@doc raw"""
    find_atoms_by_name_from_metadata(atoms::Vector{AtomMetadata}, query::Union{Vector{String}, String})::Vector{Int64}
    find_atoms_by_name_from_metadata(metadata::Metadata, query::Union{Vector{String}, String})::Vector{Int64}

Return an array of the *global* indices of all atoms whose `name` is in the `query` array.

# Examples
```julia-repl
julia> Common.find_atoms_by_name_from_metadata(metadata.atoms, ["CA", "CB"])
[5, 6, 11, 12]

julia> Common.find_atoms_by_name_from_metadata(metadata, "N")
[1, 10, 20]
```
"""
function find_atoms_by_name_from_metadata(atoms::Vector{AtomMetadata}, query::Vector{String})::Vector{Int64}

    # This function will iterate over each atom metadata, returning the *global* indices of all atoms whose name is in query.
    indices  = Vector{Int64}()
    for (index, atom) in enumerate(atoms)
        if atom.name in query
            push!(indices, index)
        end
    end
    return indices
end
find_atoms_by_name_from_metadata(atoms::Vector{AtomMetadata}, query::String)::Vector{Int64} = find_atoms_by_name_from_metadata(atoms, [query])
find_atoms_by_name_from_metadata(metadata::Metadata, query::Vector{String})::Vector{Int64} = find_atoms_by_name_from_metadata(metadata.atoms, query)
find_atoms_by_name_from_metadata(metadata::Metadata, query::String)::Vector{Int64} = find_atoms_by_name_from_metadata(metadata.atoms, [query])

@doc raw"""
    find_backbone_indices_from_metadata(atoms::Vector{AtomMetadata})::Vector{Int64}
    find_backbone_indices_from_metadata(metadata::Metadata)::Vector{Int64}

Return an array of the *global* indices of all atoms that constitute the backbone of the molecule ("N", "CA", "C").

# Examples
```julia-repl
julia> Common.find_backbone_indices_from_metadata(metadata.atoms)
[1, 5, 10, 12, 14, 21]
```
"""
function find_backbone_indices_from_metadata(atoms::Vector{AtomMetadata})::Vector{Int64}
    return find_atoms_by_name_from_metadata(atoms, ["N", "CA", "C"])
end
find_backbone_indices_from_metadata(metadata::Metadata)::Vector{Int64} = find_backbone_indices_from_metadata(metadata.atoms)

@doc raw"""
    find_ca_indices_from_metadata(atoms::Vector{AtomMetadata})::Vector{Int64}
    find_ca_indices_from_metadata(metadata::Metadata)::Vector{Int64}

Return an array of the *global* indices of all atoms that constitute the alpha carbons of the molecule ("CA").

# Examples
```julia-repl
julia> Common.find_ca_indices_from_metadata(metadata)
[10, 14, 25]
```
"""
function find_ca_indices_from_metadata(atoms::Vector{AtomMetadata})::Vector{Int64}
    return find_atoms_by_name_from_metadata(atoms, "CA")
end
find_ca_indices_from_metadata(metadata::Metadata)::Vector{Int64} = find_ca_indices_from_metadata(metadata.atoms)


@doc raw"""
    find_intra_residue_movables(atoms::Vector{AtomMetadata}, a3::Int64, current_res::Int64, current_list::Vector{Int64}, exclude::Int64)::Vector{Int64}

Return an array of the *global* indices of all atoms that constitute this dihedral `movable` list. 

# Arguments
- `atoms::Vector{AtomMetadata}`: list of all atoms metadata.
- `a3::Int64`: *global* index of the a3 atom of this dihedral.
- `current_res::Int64`: Index of the residue this dihedral belongs to.
- `current_list::Vector{Int64}`: All identified movables will be appended to this list. Used to make sure no duplicates are indentified.
- `exclude::Int64`: *Global* index. Exclude this atom when recursively searching for movables. 

# Examples
```julia-repl
julia> Common.find_intra_residue_movables(metadata.atoms, 5, 1, [5], 1)
[6, 7, 8, 9, 10, 11, 12, 13]
```
"""
function find_intra_residue_movables(atoms::Vector{AtomMetadata}, a3::Int64, current_res::Int64, current_list::Vector{Int64}, exclude::Int64)::Vector{Int64}

    # This function receives the full atoms metadata list so it can look up the current atom `res_num` and make sure the recursive search does not extend to other
    # residues AND to search for the `connects` records. The recursive search is also stopped if the connect is already in the `movables` list or is equal to the
    # exclude atom (that identifies the previously looked up atom in the chain).
    if atoms[a3].connects != nothing
        for connect in atoms[a3].connects
            (connect in current_list || atoms[connect].res_num != current_res || connect == exclude) ? continue : nothing

            push!(current_list, connect)
            current_list = find_intra_residue_movables(atoms, connect, current_res, current_list, a3)
        end
    end
    return current_list
end


@doc raw"""
    compile_backbone_dihedrals_from_metadata(atoms::Vector{AtomMetadata}[, residues::Vector{Residue}, ignore_omega::Bool = true])::Vector{Dihedral}
    compile_backbone_dihedrals_from_metadata(metadata::Metadata[, residues::Vector{Residue}, ignore_omega::Bool = true])::Vector{Dihedral}

Return a dihedrals (PHI, PSI and OMEGA) list, compiling the avaliable information from the metadata and residues list.
If `ignore_omega` flag is set to true, returns only the PHI and PSI dihedrals.

# Examples
```julia-repl
julia> bb_dihedrals = Common.compile_backbone_dihedrals_from_metadata(metadata.atoms, residues, false)
(...)

julia> bb_dihedrals = Common.compile_backbone_dihedrals_from_metadata(metadata, residues)
(...)

julia> bb_dihedrals = Common.compile_backbone_dihedrals_from_metadata(metadata)
(...)
```
See also: [`load_topology`](@ref) [`compile_sidechain_dihedrals_from_metadata`](@ref)
"""
function compile_backbone_dihedrals_from_metadata(atoms::Vector{AtomMetadata}, residues::Vector{Residue}, ignore_omega::Bool = true)::Vector{Dihedral}

    # This function iterates over the backbone atoms' indices (found using find_backbone_atoms_from_metadata() function), grouping each 4 atoms in order:
    # PHI -> OMEGA -> PSI -> (...)
    # Movables are found using find_intra_residue_movables() function and the residue number is defined as the residue in position `a3.res_num`.
    # If `ignore_omega` flag is set to true, will only return the PHI and PSI dihedrals from the backbone.
    dihedrals = Vector{Dihedral}()
    ignore_omega ? dtypes = [psi, phi] : dtypes = [psi, omega, phi]
    dtype_index::Int64 = 1
    ign::Int64 = 0
    
    indices = find_backbone_indices_from_metadata(atoms)
    for i in 1:(length(indices) - 3)
        if ignore_omega && i == 2 + 3 * ign
            ign += 1
            continue
        end
        # Define a3 for easiness
        a3 = atoms[indices[i + 2]]
        # Find this dihedral movables.
        movables = find_intra_residue_movables(atoms, indices[i + 2], a3.res_num, Int64[indices[i + 2]], indices[i + 1])
        # Add this dihedral to the dihedrals list.
        push!(dihedrals, Dihedral(indices[i], indices[i + 1], indices[i + 2], indices[i + 3], sort(movables), residues[a3.res_num], dtypes[dtype_index]))
        # If the iteration reached the end of the d_types list, return to the begining.
        dtype_index % length(dtypes) == 0 ? dtype_index = 1 : dtype_index += 1
    end
    return dihedrals
end

# Redundant functions
function compile_backbone_dihedrals_from_metadata(metadata::Metadata, residues::Vector{Residue}, ignore_omega::Bool = true)::Vector{Dihedral}
    return compile_backbone_dihedrals_from_metadata(metadata.atoms, residues, ignore_omega)
end

function compile_backbone_dihedrals_from_metadata(atoms::Vector{AtomMetadata}, ignore_omega::Bool = true)::Vector{Dihedral}
    residues = compile_residues_from_metadata(atoms)
    return compile_backbone_dihedrals_from_metadata(atoms, residues, ignore_omega)
end

function compile_backbone_dihedrals_from_metadata(metadata::Metadata, ignore_omega::Bool = true)::Vector{Dihedral}
    residues = compile_residues_from_metadata(metadata.atoms)
    return compile_backbone_dihedrals_from_metadata(metadata.atoms, residues, ignore_omega)
end


@doc raw"""
    compile_sidechain_dihedrals_from_metadata(atoms::Vector{AtomMetadata}[, residues::Vector{Residue}])::Vector{Dihedral}
    compile_sidechain_dihedrals_from_metadata(metadata::Metadata[, residues::Vector{Residue}])::Vector{Dihedral}

Return a dihedrals (CHI1, CHI2, etc) list, compiling the avaliable information from the metadata and residues list.

# Examples
```julia-repl
julia> bb_dihedrals = Common.compile_sidechain_dihedrals_from_metadata(metadata.atoms, residues)
(...)

julia> bb_dihedrals = Common.compile_sidechain_dihedrals_from_metadata(metadata, residues)
(...)

julia> bb_dihedrals = Common.compile_sidechain_dihedrals_from_metadata(metadata)
(...)
```
See also: [`load_topology`](@ref) [`compile_backbone_dihedrals_from_metadata`](@ref)
"""
function compile_sidechain_dihedrals_from_metadata(atoms::Vector{AtomMetadata}, residues::Vector{Residue})::Vector{Dihedral}

    # This function iterates over each residue (ignoring Proline) and identifies the side-chain dihedrals. The atoms list (containing each atom metadata)
    # is necessary to identify atom's names.
    dihedrals = Vector{Dihedral}()
    
    for residue in residues
        # Ignore Proline
        residue.name == "PRO" ? continue : nothing
        # Create a new atom_name -> atom_index conversion dictionary for each residue iterated over.
        residue_atoms = Dict{String, Int64}()
        for atom_index in residue.atoms
            # Don't add hidrogens to the conversion dictionary since they will not be accounted for when identifying dihedrals.
            startswith(atoms[atom_index].name, "H") ? continue : nothing
            # The conversion dictionary basically gives the atom *global* index based on it's name.
            residue_atoms[atoms[atom_index].name] = atom_index
        end

        # Start the dihedral path
        path = ["N", "CA"]
        # Based on the `residue.name`, the possible path taken by the iteration is different. This accounts for cyclic aminoacids.
        possible_queries = residue.name in ("HIE", "PHE", "TRP", "TYR") ? "BGD" : "BGDEZH"
        # Join all atom names in this residue in a single string, dividing atom names with ":". This will allow identification of the dihedral using regular expressions.
        atnames = string(":", join(map(at -> atoms[at].name, residue.atoms), ":"), ":")
        # Iterate over the given possible path for the side-chain.
        for (chi_number, query) in enumerate(possible_queries)
            # The regular expression tries to identify the current query in the list of atom names. This is an overpowered version of looking for a key in the dictionary
            # since it account for possible numbers in the atom name and different elements (C vs N for example).
            m = match(Regex("(.*):(?P<name>[A-GI-Z]$query[1]?):(.*)"), atnames)
            # If the current query is found, add this atom name to the dihedral path being built.
            if m != nothing
                push!(path, (m[:name]))
            end
            # If the dihedral path reaches 4 atoms (meaning a dihedral was identified), add a new Dihedral to the dihedrals list.
            if length(path) == 4
                # The dihedral type is defined in order, starting in chi1, chi2, etc. Since the `enumerate` function, in Julia, starts in 1 and the first dihedral will
                # only be identified in the second iteration of the loop (since on the first iteration only 3 atoms of the path have been identified), `$(chi_number - 1)`
                # is used.
                dtype = eval(Symbol("chi$(chi_number - 1)"))
                # The previously defined atom_name -> atom_index conversion dictionary is employed to define the dihedral by its *global* indices.
                atom_indices = map(at -> residue_atoms[at], path)
                # a3 is defined for easiness.
                a3 = atom_indices[3]
                # The dihedral movables are identified recursively using find_intra_residue_movables() function.
                movables = find_intra_residue_movables(atoms, a3, atoms[a3].res_num, Int64[a3], atom_indices[2])
                push!(dihedrals, Dihedral(atom_indices[1], atom_indices[2], a3, atom_indices[4], sort(movables), residue, dtype))
                # Delete the begining of the dihedral path. If a new atom is identified in the side-chain, the a2 of this dihedral becomes a1, etc
                deleteat!(path, 1)
            end
        end
    end
    return dihedrals
end

# Redundant functions
function compile_sidechain_dihedrals_from_metadata(metadata::Metadata, residues::Vector{Residue})::Vector{Dihedral}
    return compile_sidechain_dihedrals_from_metadata(metadata.atoms, residues)
end

function compile_sidechain_dihedrals_from_metadata(atoms::Vector{AtomMetadata})::Vector{Dihedral}
    residues = compile_residues_from_metadata(atoms)
    return compile_sidechain_dihedrals_from_metadata(atoms, residues)
end

function compile_sidechain_dihedrals_from_metadata(metadata::Metadata)::Vector{Dihedral}
    residues = compile_residues_from_metadata(metadata.atoms)
    return compile_sidechain_dihedrals_from_metadata(metadata.atoms, residues)
end