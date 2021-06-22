"""
    load([::Type{T}], filename::AbstractString; bonds_by_distance::Bool = false) where {T <: AbstractFloat}

Load the given `filename` into a [`Pose`](@ref), parametrized by `T`. If this is
not provided, the default `ProtoSyn.Units.defaultFloat` is used instead. The
file format is infered from the extension (Supported: .pdb, .yml). If
`bonds_by_distance` is set to `true` (`false`, by default), the CONECT records
will be complemented with bonds infered by distance. The threshold distances for
each pair of atoms is defined in `ProtoSyn.bond_lengths`. Infers parenthood and
ascedence from bonds (N-[`Residue`](@ref) instances have the connected
C-[`Residue`](@ref) as child).

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
function load(::Type{T}, filename::AbstractString; bonds_by_distance::Bool = false) where {T <: AbstractFloat}

    if !bonds_by_distance && ProtoSyn.verbose.mode
        @info "Flag `bonds_by_distance` is set to False. Make sure the loaded $filename file has connect records."
    end

    pose = ProtoSyn.load(T, filename, bonds_by_distance = bonds_by_distance)

    # Set parenthood of residues
    for segment in eachsegment(pose.graph)

        n_residues = ProtoSyn.count_residues(segment)
        for residue_index in 2:n_residues
            residue = segment[residue_index]
            popparent!(residue) # Was root before.
            C_bond_index = findfirst(x -> x.symbol == "C", residue["N"].bonds)
            C_bond_index == nothing && error("No atom \"C\" found in atom $(residue["N"]) bonds list ($(residue["N"].bonds)).")
            parent_residue = residue["N"].bonds[C_bond_index].container
            setparent!(residue, parent_residue)
        end
    end

    return pose
end

load(filename::AbstractString; bonds_by_distance::Bool = false) = begin
    Peptides.load(ProtoSyn.Units.defaultFloat, filename, bonds_by_distance = bonds_by_distance)
end