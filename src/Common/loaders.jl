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
    xyz   = Vector{Array{Float64, 2}}()
    atoms = Vector{AtomMetadata}()

    #Read file (from file name) and recover XYZ and ATOM NAMES
    open(i_file, "r") do f
        for (index, line) in enumerate(eachline(f))
            elem = split(line)
            if length(elem) > 3 && index > 2
                push!(xyz, map(x -> parse(Float64, x), [elem[4] elem[5] elem[6]]))
                push!(atoms, AtomMetadata(string(elem[2]), elem=string(elem[2]), res_num=parse(Int64, strip(line[1:5])), res_name=string(strip(line[6:8]))))
            end
        end
    end

    n = length(xyz)
    return Common.State(n, Common.Energy(), vcat(xyz...), zeros(n, 3)), Metadata(atoms = atoms)
end


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
function load_from_pdb(i_file::String, compile_metadata = true)

    xyz   = Vector{Array{Float64, 2}}()
    if compile_metadata atoms = Vector{AtomMetadata}() end

    open(i_file, "r") do f
        for line in eachline(f)
            if startswith(line, "ATOM")
                push!(xyz, map(x -> 0.1*parse(Float64, x), [line[31:38] line[39:46] line[47:54]]))
                if compile_metadata
                    push!(atoms, AtomMetadata(
                        index    = parse(Int64, strip(line[6:11])),
                        name     = string(strip(line[14:16])),
                        elem     = string(strip(line[77:78])),
                        res_num  = parse(Int64, line[23:26]),
                        res_name = string(strip(line[18:20])),
                        chain_id = string(line[22]),
                        connects = nothing))
                end
            elseif compile_metadata && startswith(line, "CONECT")
                elem = split(line)
                atoms[parse(Int64, elem[2])].connects = map(x -> parse(Int64, x), elem[3:end])
            end
        end
    end
    n = length(xyz)

    if compile_metadata
        metadata = Metadata(atoms = atoms)
        compile_residue_metadata!(metadata)
        compile_dihedral_metadata!(metadata)
    else
        metadata = Metadata()
    end

    return State(n, Energy(), vcat(xyz...), zeros(n, 3)), metadata
end


@doc raw"""
    load_metadata_from_json(json_file::String)

Parse a JSON file containing the dihedral and residue topology. Return a [`Dihedral`](@ref) array
and a [`Residue`](@ref) array.

# Examples
```julia-repl
julia> Mutators.Diehdral.load_topology("mc.json")
(ProtoSyn.Mutators.Dihedral.NewDihedral[...], ProtoSyn.Common.Residue[...])
```
See also: [`Aux.read_JSON`](@ref)
"""
function load_metadata_from_json(json_file::String)

    p::Dict{String, Any} = Aux.read_JSON(json_file)
    tmp_residues = Dict(d["n"] => Residue(convert(Vector{Int64}, d["atoms"]), d["next"], d["type"]) for d in p["residues"])
    
    str2enum = Dict(string(s) => s for s in instances(DIHEDRAL.TYPE))
    
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
See also: [`load_metadata_from_json`](@ref)
"""
function compile_residue_metadata!(metadata::Metadata)

    # Verify input
    if length(metadata.atoms) <= 0 error("Metadata needs to have AtomMetadata information in order to compile residues.") end

    # This function will iterate over each atom metadata in reverse order. This allows to define the `residue.next` as the previously defined residue in the iteration,
    # since it is being looked up in reverse. The first residue defined does not have a `next`, so it is exeptionally set to `nothing`. A new residue is identified
    # whenever the atom being iterated over has a different `res_num` than the previous one. Newly identified residues are appended to the begining of the list, so
    # they are returned in order.
    residue_metadata = Vector{Residue}()
    
    # Iterate over the atom list in reverse order
    residues = Dict{Int64, Vector{AtomMetadata}}()
    for atom in metadata.atoms
        if !(atom.res_num in keys(residues))
            residues[atom.res_num] = Vector{AtomMetadata}()
        end
        push!(residues[atom.res_num], atom)
    end

    # Iterate over the reversed residues list
    for residue in reverse(sort(collect(values(residues)), by = atoms -> atoms[1].index))
        # Exceptionally, the "first" residue (actually last, since it's in reverse order) will have `nothing` as its `next`.
        if length(residue_metadata) == 0
            insert!(residue_metadata, 1, Residue([atom.index for atom in residue], nothing, residue[1].res_name, SS.COIL))
            for atom in residue atom.residue = residue_metadata[1] end
        else
            # Add this residue to the residues list, appending it to the begining of the list.
            # All other residues have the previously defined residue as its `next`.
            insert!(residue_metadata, 1, Residue([atom.index for atom in residue], residue_metadata[1], residue[1].res_name, SS.COIL))
            for atom in residue atom.residue = residue_metadata[1] end
        end
    end
    # Append the first residue to the begining of the list
    metadata.residues = residue_metadata
end


@doc raw"""
    compile_dihedrals_from_metadata(atoms::Vector{AtomMetadata})::Vector{Dihedral}
    compile_dihedrals_from_metadata(metadata::Metadata)::Vector{Dihedral}

Return a dihedrals list, compiling the avaliable information from the metadata.

# Examples
```julia-repl
julia> dihedrals = Common.compile_dihedrals_from_metadata(metadata.atoms)
(...)

julia> dihedrals = Common.compile_dihedrals_from_metadata(metadata)
(...)
```
See also: [`load_metadata_from_json`](@ref)
"""
function compile_dihedral_metadata!(metadata::Metadata)

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

    #Verify input
    if length(metadata.atoms) <= 0 error("Metadata needs to have AtomMetadata information in order to compile dihedrals.") end

    # This function iterates over each residue (ignoring Proline) and identifies the side-chain dihedrals. The atoms list (containing each atom metadata)
    # is necessary to identify atom's names.
    dihedrals = Vector{Dihedral}()
    
    # Identify residues from atom data
    residues = Dict{Int64, Vector{AtomMetadata}}()
    for atom in metadata.atoms
        if !(atom.res_num in keys(residues))
            residues[atom.res_num] = Vector{AtomMetadata}()
        end
        push!(residues[atom.res_num], atom)
    end

    for residue in sort(collect(values(residues)), by = atoms -> atoms[1].index)
        # Ignore Proline
        if residue[1].res_name == "PRO" continue end
        # Create a new atom_name -> atom_index conversion dictionary for each residue iterated over.
        name2index = Dict{String, Int64}()
        for atom in residue
            # Don't add hidrogens to the conversion dictionary since they will not be accounted for when identifying dihedrals.
            if startswith(atom.name, "H") continue end
            # The conversion dictionary basically gives the atom *global* index based on it's name.
            name2index[atom.name] = atom.index
        end

        # Identify Backbone Dihedrals
        n, ca, c = filter(atom -> atom.name in ["N", "CA", "C"], residue)
        # PHI
        prev_c = filter(connect -> metadata.atoms[connect].name == "C", n.connects)
        if length(prev_c) > 0
            movables = find_intra_residue_movables(metadata.atoms, ca.index, ca.res_num, Int64[ca.index], n.index)
            push!(dihedrals, Dihedral(prev_c[1], n.index, ca.index, c.index, sort(movables), ca.residue, DIHEDRAL.phi))
        end
        # PSI
        next_n = filter(connect -> metadata.atoms[connect].name == "N", c.connects)
        if length(next_n) > 0
            movables = find_intra_residue_movables(metadata.atoms, c.index, c.res_num, Int64[c.index], ca.index)
            push!(dihedrals, Dihedral(n.index, ca.index, c.index, next_n[1], sort(movables), c.residue, DIHEDRAL.psi))
            # OMEGA
            next_n = metadata.atoms[next_n[1]]
            next_ca = filter(connect -> metadata.atoms[connect].name == "CA", next_n.connects)
            if length(next_ca) > 0
                movables = find_intra_residue_movables(metadata.atoms, next_n.index, next_n.res_num, Int64[next_n.index], c.index)
                push!(dihedrals, Dihedral(ca.index, c.index, next_n.index, next_ca[1], sort(movables), next_n.residue, DIHEDRAL.omega))
            end
        end

        # Identify Side-chain Dihedrals
        # Start the dihedral path
        path = ["N", "CA"]
        # Based on the `residue.name`, the possible path taken by the iteration is different. This accounts for cyclic aminoacids.
        possible_queries = residue[1].res_name in ("HIE", "PHE", "TRP", "TYR") ? "BGD" : "BGDEZH"
        # Join all atom names in this residue in a single string, dividing atom names with ":". This will allow identification of the dihedral using regular expressions.
        atnames = string(":", join(map(atom -> atom.name, residue), ":"), ":")
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
                # The previously defined atom_name -> atom_index conversion dictionary is employed to define the dihedral by its *global* indices.
                a1, a2, a3, a4 = map(atom_name -> name2index[atom_name], path)
                # The dihedral movables are identified recursively using find_intra_residue_movables() function.
                movables = find_intra_residue_movables(metadata.atoms, a3, metadata.atoms[a3].res_num, Int64[a3], a2)
                # The dihedral type is defined in order, starting in chi1, chi2, etc. Since the `enumerate` function, in Julia, starts in 1 and the first dihedral will
                # only be identified in the second iteration of the loop (since on the first iteration only 3 atoms of the path have been identified), `chi_number - 1`
                # is used.
                push!(dihedrals, Dihedral(a1, a2, a3, a4, sort(movables), residue[1].residue, DIHEDRAL.TYPE(chi_number - 1)))
                # Delete the begining of the dihedral path. If a new atom is identified in the side-chain, the a2 of this dihedral becomes a1, etc
                deleteat!(path, 1)
            end
        end
    end
    metadata.dihedrals = dihedrals
end