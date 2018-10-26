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
function load_from_pdb(i_file::String)::State

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
    return State(n, Energy(), vcat(xyz...), zeros(n, 3), Metadata(metadata))
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

    # println(p["tmp_tmp_residues"][1]["atoms"], typeof(p["tmp_residues"][1]["atoms"]))
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