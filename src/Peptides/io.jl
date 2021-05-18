"""
    load([::Type{T}], filename::AbstractString; bonds_by_distance::Bool = false, verbose::Bool = true) where {T <: AbstractFloat}

Load the given `filename` into a pose, parametrized by T. If this is not
provided, the default ProtoSyn.Units.defaultFloat is used. The file format is
infered from the extension (Supported: .pdb, .yml). If `bonds_by_distance` is
set to true (false, by default), the CONECT records will be complemented with
bonds infered by distance. The threshold distances for each pair of atoms is
defined in ProtoSyn.bond_lengths. Infer parenthood and ascedence from bonds. If
`verbose` is set to `true` (is, by default), print the loading status.

!!! note
    This function is an overload of `ProtoSyn.load`.

# Examples
```
julia> ProtoSyn.Peptides.load("1ctf.pdb", bonds_by_distance = true)
...
```
"""
function load(::Type{T}, filename::AbstractString; bonds_by_distance::Bool = false, verbose::Bool = true) where {T <: AbstractFloat}

    if !bonds_by_distance && verbose
        @info "Flag `bonds_by_distance` is set to False. Make sure the loaded $filename file has connect records."
    end

    pose = ProtoSyn.load(T, filename, bonds_by_distance = bonds_by_distance)

    # Set parenthood pf residues
    for segment in eachsegment(pose.graph)
        setparent!(segment[1][1], ProtoSyn.root(pose.graph))

        n_residues = ProtoSyn.count_residues(segment)
        for residue_index in 2:n_residues
            residue = segment[residue_index]
            popparent!(residue)
            C_bond_index = findfirst(x -> x.symbol == "C", residue["N"].bonds)
            parent_residue = residue["N"].bonds[C_bond_index].container
            setparent!(residue, parent_residue)
        end
    end

    # Set parenthood of atoms
    atoms   = collect(eachatom(pose.graph))
    n_atoms = length(atoms)
    visited = ProtoSyn.Mask{Atom}(n_atoms)
    for (i, atom_i) in enumerate(atoms)
        visited[atom_i.index] = true
        for atom_j in atom_i.bonds
            if visited[atom_j.index]
            else
                !hasparent(atom_j) && ProtoSyn.setparent!(atom_j, atom_i)
            end
        end
    end

    reindex(pose.graph)
    ProtoSyn.request_c2i!(pose.state)
    sync!(pose)

    return pose
end

load(filename::AbstractString; bonds_by_distance::Bool = false) = begin
    Peptides.load(ProtoSyn.Units.defaultFloat, filename, bonds_by_distance = bonds_by_distance)
end