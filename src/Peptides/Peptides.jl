@doc """
    Peptides

The Peptides modules introduces `Calculators`, `Mutators` and `Drivers` (among
other methods) specific for proteins and peptides.
"""
module Peptides

using ..ProtoSyn
using ..ProtoSyn.Builder

# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname)
end

include("constants.jl")
include("Calculators/Calculators.jl")
include("Rotamers/Rotamers.jl")

export grammar

"""
    grammar([::Type{T}]) where {T <: AbstractFloat}

Build a `LGrammar` for peptides, taking as variables the fragments in the
default resource directory. If the option type T is not provided, the default
ProtoSyn float value will be used. The returned LGrammar is required for
building peptides from fragments.

# Examples
```julia-repl
julia> g = Peptides.grammar();
julia> pose = Builder.build(grammar, seq"AAGASTASSE")
...
```
"""
function grammar(::Type{T}) where {T <: AbstractFloat}
    filename = joinpath(resource_dir, "grammars.yml")
    Builder.fromfile(T, filename, "peptide")
end

grammar() = grammar(ProtoSyn.Units.defaultFloat)


"""
"""
function build(grammar::LGrammar{T}, derivation, ss::NTuple{3,Number} = SecondaryStructure[:linear]) where {T <: AbstractFloat}

    pose = Builder.build(grammar, derivation)
    Peptides.setss!(pose, ss)
    sync!(pose)
    pose
end


"""
"""
function load(::Type{T}, filename::AbstractString; bonds_by_distance::Bool = false) where {T <: AbstractFloat}

    pose = ProtoSyn.load(T, filename)
    
    for segment in eachsegment(pose.graph)
        setparent!(segment[1][1], ProtoSyn.origin(pose.graph))
        
        n_residues = ProtoSyn.count_residues(segment)
        for residue_index in 2:n_residues
            popparent!(segment[residue_index])
            setparent!(segment[residue_index], segment[residue_index - 1])
        end
    end

    if bonds_by_distance
        dm        = ProtoSyn.Calculators.full_distance_matrix(pose)
        threshold = T(0.1)
    end

    atoms   = collect(eachatom(pose.graph))
    n_atoms = length(atoms)
    visited = ProtoSyn.Mask{Atom}(n_atoms)
    for (i, atom_i) in enumerate(atoms)
        visited[atom_i.index] = true
        if bonds_by_distance
            for (j, atom_j) in enumerate(atoms)
                i == j && continue
                atom_j = atoms[j]
                atom_j in atom_i.bonds && continue
                putative_bond = "$(atom_i.symbol)$(atom_j.symbol)"
                !(putative_bond in keys(Peptides.bond_lengths)) && continue
                d = Peptides.bond_lengths[putative_bond]
                d += d * threshold
                dm[i, j] < d && ProtoSyn.bond(atom_i, atom_j)
            end
        end
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

    pose
end

load(filename::AbstractString; bonds_by_distance::Bool = false) = begin
    Peptides.load(ProtoSyn.Units.defaultFloat, filename, bonds_by_distance = bonds_by_distance)
end




function uncap!(pose::Pose, residue::Residue)

    @assert length(residue["N"].children) > 2 "Residue $residue doesn't seem to be capped."

    for atom in ProtoSyn.gather(an"H\d+"r(pose), residue)
        pop!(pose, atom)
    end

    # b = 0.9795641888105121
    # θ = 2.095336203538727
    # ϕ = 0.0001555866587524742

    H = Atom!(residue, "H", -3, -1, "H")
    setparent!(H, residue["N"])
    ProtoSyn.bond(residue["N"], H)

    reindex(pose.graph)

    insert!(pose.state, H.index, State(1))
    pose.state[H].b = 0.9795641888105121
    pose.state[H].θ = deg2rad(120)
    pose.state[H].ϕ = deg2rad(180)

    ProtoSyn.request_i2c(pose.state)
end


"""
    sequence(container::ProtoSyn.AbstractContainer)::String

Return the sequence of aminoacids (in 1 letter mode) of the given container/pose
as a string.

# Examples
```julia-repl
julia> sequence(pose)
"AAGASTASSE"
```
"""
function sequence(container::ProtoSyn.AbstractContainer)::String

    sequence = ""
    for residue in eachresidue(container)
        try
            sequence *= three_2_one[residue.name]
        catch KeyError
            sequence *= '?'
        end
    end

    return sequence
end

sequence(pose::Pose) = sequence(pose.graph)


include("methods.jl")
include("Mutators/Mutators.jl")
include("Drivers/drivers.jl")

end