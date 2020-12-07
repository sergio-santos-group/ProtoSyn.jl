"""
    load([::Type{T}], filename::AbstractString; bonds_by_distance::Bool = false) where {T <: AbstractFloat}

Load the given `filename` into a pose, parametrized by T. If this is not provided,
the default ProtoSyn.Units.defaultFloat is used. The file format is infered from
the extension (Supported: .pdb, .yml). If `bonds_by_distance` is set to true
(false, by default), the CONECT records will be complemented with bonds infered
by distance. The threshold distances for each pair of atoms is defined in
ProtoSyn.bond_lengths. Infer parenthood and ascedence from bonds.

# Examples
```julia-repl
julia> Peptides.load("1ctf.pdb", bonds_by_distance = true)
...
```
"""
function load(::Type{T}, filename::AbstractString; bonds_by_distance::Bool = false) where {T <: AbstractFloat}

    pose = ProtoSyn.load(T, filename, bonds_by_distance = bonds_by_distance)

    for segment in eachsegment(pose.graph)
        setparent!(segment[1][1], ProtoSyn.origin(pose.graph))

        n_residues = ProtoSyn.count_residues(segment)
        for residue_index in 2:n_residues
            popparent!(segment[residue_index])
            C_bond_index = findfirst(x -> x.symbol == "C", residue["N"].bonds)
            parent_residue = residue["N"].bonds[C_bond_index].container
            setparent!(segment[residue_index], parent_residue)
        end
    end

    atoms   = collect(eachatom(pose.graph))
    n_atoms = length(atoms)
    visited = ProtoSyn.Mask{Atom}(n_atoms)
    for (i, atom_i) in enumerate(atoms)
        visited[atom_i.index] = true
        for atom_j in atom_i.bonds
            if visited[atom_j.index]
                atom_i.parent = atom_j
            else
                push!(atom_i.children, atom_j)
            end
        end
    end

    reindex(pose.graph)
    ProtoSyn.request_c2i(pose.state)
    sync!(pose)

    return pose
end

load(filename::AbstractString; bonds_by_distance::Bool = false) = begin
    Peptides.load(ProtoSyn.Units.defaultFloat, filename, bonds_by_distance = bonds_by_distance)
end