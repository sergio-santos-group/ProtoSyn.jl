"""
    load([::Type{T}], filename::AbstractString; bonds_by_distance::Bool = false; alternative_location::String = "A") where {T <: AbstractFloat}

Load the given `filename` into a [`Pose`](@ref), parametrized by `T`. If this is
not provided, the default `ProtoSyn.Units.defaultFloat` is used instead. The
file format is infered from the extension (Supported: .pdb, .yml). If
`bonds_by_distance` is set to `true` (`false`, by default), the CONECT records
will be complemented with bonds infered by distance. The threshold distances for
each pair of atoms is defined in `ProtoSyn.bond_lengths`. Infers parenthood and
ascedence from bonds (N-[`Residue`](@ref) instances have the connected
C-[`Residue`](@ref) as child). By default, and when available, ProtoSyn will use
`alternative_location` `A`, unless specified in the flag `alternative_location`.

!!! ukw "Note:"
    This function is an extension of [`ProtoSyn.load`](@ref).

# Examples
```
julia> ProtoSyn.Peptides.load("1ctf.pdb", bonds_by_distance = true)
Pose{Topology}(Topology{/2a3d:61708}, State{Float64}:
 Size: 1140
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function load(::Type{T}, filename::AbstractString; bonds_by_distance::Bool = false, alternative_location::String = "A") where {T <: AbstractFloat}

    if !bonds_by_distance && ProtoSyn.verbose.mode
        @info "Flag `bonds_by_distance` is set to False. Make sure the loaded $filename file has connect records."
    end

    pose = ProtoSyn.load(T, filename, bonds_by_distance = bonds_by_distance, alternative_location = alternative_location)

    available_aminoacids = [ProtoSyn.Peptides.one_2_three[x] for x in keys(ProtoSyn.Peptides.available_aminoacids)]
    
    # Set parenthood of residues
    for segment in eachsegment(pose.graph)

        n_residues = ProtoSyn.count_residues(segment)
        for residue_index in 2:n_residues
            residue = segment[residue_index]
            if !(residue.name in available_aminoacids)
                @warn "Found a possible ligand at residue $residue. ProtoSyn will skip it when setting up residue parenthoods."
                continue
            end
            popparent!(residue) # Was root before.
            C_bond_index = findfirst(x -> x.symbol == "C", residue["N"].bonds)
            C_bond_index == nothing && error("No atom \"C\" found in atom $(residue["N"]) bonds list ($(residue["N"].bonds)).")
            parent_residue = residue["N"].bonds[C_bond_index].container
            setparent!(residue, parent_residue)
        end
    end

    return pose
end

load(filename::AbstractString; bonds_by_distance::Bool = false, alternative_location::String = "A") = begin
    Peptides.load(ProtoSyn.Units.defaultFloat, filename, bonds_by_distance = bonds_by_distance, alternative_location = alternative_location)
end


"""
    ProtoSyn.download([::T], pdb_code::String) where {T <: AbstractFloat}

Download the PDB file (for the given PDB code) from the RCSB
Protein Data Bank into a [`Pose`](@ref). The downloaded file can be found in the current working
directory. If `T` is specified, the downloaded file will be loaded into a
[`Pose`](@ref) parametrized by `T`, otherwise uses the default
`ProtoSyn.Units.defaultFloat`. Uses the specific `Peptides.load` method.

# See also
[`load`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.download("2A3D")
```
"""
function ProtoSyn.Peptides.download(::Type{T}, pdb_code::String; bonds_by_distance::Bool = false) where {T <: AbstractFloat}
    if endswith(pdb_code, ".pdb"); pdb_code = pdb_code[1:(end - 4)]; end
    filename = pdb_code * ".pdb"
    url = "https://files.rcsb.org/download/" * filename
    download(url, filename)
    return ProtoSyn.Peptides.load(T, filename, bonds_by_distance = bonds_by_distance)
end

ProtoSyn.Peptides.download(pdb_code::String; bonds_by_distance::Bool = false) = begin
    ProtoSyn.Peptides.download(ProtoSyn.Units.defaultFloat, pdb_code, bonds_by_distance = bonds_by_distance)
end