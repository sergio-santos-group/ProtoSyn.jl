# ----------------------------------------------------------------------------------------------------------
#                                                 LOADERS

#TODO: Verify function
@doc raw"""
    load_from_gro(i_file::String)::Common.State

Return a new [`Common.State`](@ref) by loading the atom positions and names from the input .gro file.
As a default, `state.energy` is [`NullEnergy`](@ref) and `state.forces` are set to zero.

# Examples
```julia-repl
julia> Common.load_from_gro("molecule.gro")
Common.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], atnames=["C", "O"])
```
See also: [`load_from_pdb`](@ref)
"""
function load_from_gro(i_file::String)::Common.State

    #Initialize empty arrays
    xyz     = Array{Array{Float64, 2}, 1}()
    atnames = Array{String, 1}()

    #Read file (from file name) and recover XYZ and ATOM NAMES
    open(i_file, "r") do f
        for (index, line) in enumerate(eachline(f))
            elem = split(line)
            if length(elem) > 3 && index > 2
                push!(xyz, map(x -> parse(Float64, x), [elem[4] elem[5] elem[6]]))
                push!(atnames, elem[2])
            end
        end
    end

    n = length(xyz)
    return Common.State(n, Common.NullEnergy(), vcat(xyz...), zeros(n, 3), atnames)
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
function load_from_pdb(i_file::String)::State

    xyz      = Array{Array{Float64, 2}, 1}()
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
    return State(n, NullEnergy(), vcat(xyz...), zeros(n, 3), Metadata(metadata))
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

    residues = Dict(d["n"] => Residue(d["atoms"], d["next"], d["type"]) for d in p["residues"])
    
    str2enum = Dict(string(s) => s for s in instances(DIHEDRALTYPE))
    
    dihedrals = [
        Dihedral(d["a1"], d["a2"], d["a3"], d["a4"],
            d["movable"], residues[d["parent"]], str2enum[lowercase(d["type"])])
        for d in p["dihedrals"]
    ]
    
    # Set correct references for dihedrals previous and next
    for residue in values(residues)
        residue.next = get(residues, residue.next, nothing)
    end

    return dihedrals, residues
end