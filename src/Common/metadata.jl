# ----------------------------------------------------------------------------------------------------------
#                                                METADATA

@doc raw"""
    AtomMetadata(name::String[, elem::String = name, res_num::Int64 = 1, res_name::String = "UNK", chain_id::Union{String, Nothing} = nothing, connects::Union{Vector{Int64}, Nothing} = nothing])

Define an atom metadata, containing extra information pertaining the [`State`](@ref).

# Arguments
- `index::Int64`: *global* index of this atom.
- `name::String`: Name of the atom.
- `elem::String`: (Optional) Element of the atom (Default: `name`).
- `res_num::Int64`: (Optional) Number of the residue this atom belongs to (Default: 1).
- `res_name::Union{String, Nothing}`: (Optional) Name of the residue this atom belongs to (Default: "UNK").
- `residue::Union{Residue, Nothing}`: (Optional) Reference to the residue this atom belong to (Deffault: nothing).
- `chain_id::String`: (Optional) Name of the chain that contains the residue this atom belongs to (Default: nothing).
- `connects::Union{Vector{Int64}, Nothing}`: (Optional) List of *global* atom indices that this atom is connected to (Default: nothing). 

# Examples
```julia-repl
julia> AtomMetadata("H1", "H", 2, "VAL", Residue(...), "A", [4])
AtomMetadata(name=H1, elem=H, res_num=2, res_name=VAL, chain_id=A, connects=[4])

julia> AtomMetadata("H1")
AtomMetadata(name=H1, elem=H1, res_num=1, res_name=UNK, chain_id=nothing, connects=nothing)
```
"""
Base.@kwdef mutable struct AtomMetadata

    # Parameter                             Default
    index::Int64                            = -1
    name::String                            = "X"
    symbol::String                          = "X"
    residue::Union{Residue, Nothing}        = nothing
    connects::Union{Vector{Int64}, Nothing} = nothing

    # AtomMetadata(;index::Int64 = 0, name::String = "_", elem::String = name, res_num::Int64 = 1, res_name::String = "UNK", residue = nothing, chain_id::Union{String, Nothing} = nothing, connects::Union{Vector{Int64}, Nothing} = nothing) = new(index, name, elem, res_num, res_name, residue, chain_id, connects)
end
# function Base.show(io::IO, b::AtomMetadata) #SHOULD BE IMPROVED
#     if b.chain_id != nothing && b.connects != nothing
#         print(io, "AtomMetadata(index = $(b.index), name=$(b.name), elem=$(b.elem), res_num=$(b.res_num), res_name=$(b.res_name), chain_id=$(b.chain_id), connects=$(b.connects))")
#     elseif b.chain_id != nothing
#         print(io, "AtomMetadata(index = $(b.index), name=$(b.name), elem=$(b.elem), res_num=$(b.res_num), res_name=$(b.res_name), chain_id=$(b.chain_id), connects=nothing)")
#     elseif b.connects != nothing
#         print(io, "AtomMetadata(index = $(b.index), name=$(b.name), elem=$(b.elem), res_num=$(b.res_num), res_name=$(b.res_name), chain_id=nothing, connects=$(b.connects))")
#     else
#         print(io, "AtomMetadata(index = $(b.index), name=$(b.name), elem=$(b.elem), res_num=$(b.res_num), res_name=$(b.res_name), chain_id=nothing, connects=nothing)")
#     end
# end

function Base.show(io::IO, b::AtomMetadata)
    print(io, string(typeof(b)))
    for p in fieldnames(typeof(b))
        v = getproperty(b,p)
        s = v != nothing ? v : "nothing"
        print(io, "\n   $(String(p)) = $v")
    end
end

function iterate_by_residue(atoms::Vector{AtomMetadata})::Vector{Vector{AtomMetadata}}

    residue::Vector{AtomMetadata} = AtomMetadata[]
    chain::Vector{Vector{AtomMetadata}} = []
    old_res_num = 1
    for atom in atoms
        if atom.res_num != old_res_num
            old_res_num = atom.res_num
            push!(chain, residue)
            residue = AtomMetadata[]
        end
        push!(residue, atom)
    end
    push!(chain, residue)
    return chain
end


@doc raw"""
    SecondaryStructureMetadata(ss_type::SS.TYPE, name::String, i_res_name::String, i_res_num::Int64, f_res_name::String, f_res_num::Int64, conf::Int64)

Define a secondary structure metadata, containing extra information pertaining the [`State`](@ref).

# Arguments
- `ss_type::SS.TYPE`: Type of the SecondaryStructure (SS.HELIX, SS.SHEET, ...)
- `name::String`: Name of the structure.
- `i_res_name::String`: Name of the starting residue.
- `i_res_num::Int64`: Index of the starting residue.
- `f_res_name::String`: Name of the final residue.
- `f_res_num::Int64`: Index of the final residue.
- `conf::Int64`: Conformation of the secondary structure. See PDB FORMAT standards.

# Examples
```julia-repl
julia> SecondaryStructureMetadata(SS.HELIX, "HA", "V", 4, "A", 7, 1)
SecondaryStructureMetadata(ss_type=HELIX, name=HA, V-4 <-> A-7, conf=1)
```
"""
mutable struct SecondaryStructureMetadata

    ss_type::SS.TYPE
    name::String
    i_res_name::String
    i_res_num::Int64
    f_res_name::String
    f_res_num::Int64
    conf::Int64

end
Base.show(io::IO, b::SecondaryStructureMetadata) = print(io, "SecondaryStructureMetadata(ss_type=$(b.ss_type), name=$(b.name), $(b.i_res_name)-$(b.i_res_num) <-> $(b.f_res_name)-$(b.f_res_num), conf=$(b.conf))")


@doc raw"""
    BlockMetadata(atoms::Vector{Int64}, pivot::Int64, range_left::Float64, connector_left::Int64, connector_right::Int64)

Define a block, containing all the necessary information for [`BlockrotMutator`] (@ref Blockrot) movements.

# Arguments
- `atoms::Vector{Int64}`: List of *global* indexes of all atoms contained in the block.
- `pivot::Int64`: *Global* index of the center atom of the block.
- `range_left::Float64`: Range (in nm), of the left side coil.
- `connector_left::Int64`: *Global* index of the left side connector of the block (N).
- `connector_right::Int64`: *Global* index of the right side connector of the block (C).

# Examples
```julia-repl
julia> Common.BlockMetadata([1, 2, 3], 2, 0.8, 1, 3)
BlockMetadata(atoms=1<->3, pivot=2, range_left=0.8, connector_left=1, connector_right=3)
```
See also: 
"""
mutable struct BlockMetadata

    atoms::Vector{Int64}
    pivot::Int64
    range_left::Float64
    connector_left::Int64
    connector_right::Int64
end
Base.show(io::IO, b::BlockMetadata) = print(io, "BlockMetadata(atoms=$(b.atoms[1])<->$(b.atoms[length(b.atoms)]), pivot=$(b.pivot), range_left=$(b.range_left), connector_left=$(b.connector_left), connector_right=$(b.connector_right))")


struct SpanningTree
    size::Int
    table::Vector{Vector{Int}}
    stack::Vector{Int}
    indices::Vector{Int}
    visited::BitVector
end
SpanningTree(table::Vector{Vector{Int}}) = begin
    n = length(table)
    visited = BitVector(undef, n)
    SpanningTree(n, table, zeros(Int, n), zeros(Int, n), visited)
end

function getindex(tree::SpanningTree, pivot::Int)
    if (pivot < 1) || (pivot > tree.size)
        return 0
    end
    tree.visited .= false

    j = 1
    ptr = 1
    tree.stack[ptr] = pivot
    @inbounds while ptr > 0
        id = tree.stack[ptr]
        ptr -= 1
        if !tree.visited[id]
            tree.visited[id] = true
            tree.indices[j] = id
            j += 1
            for i in tree.table[id]
                ptr += 1
                tree.stack[ptr] = i
            end
        end
    end
    j-1
end

# @doc raw"""
#     Metadata([, atoms::Vector{AtomMetadata} = [], ss::Vector{SecondaryStructureMetadata} = [], residues::Vector{Residue} = [], dihedrals::Vector{Dihedral} = [], blocks::Vector{BlockMetadata} = [], sidechains::Vector{SidechainMetadata} = []])

# Define the state metadata, containing extra information regarding the atoms and secondary structure of the system.

# # Examples
# ```julia-repl
# julia> Metadata(atoms, ss, residues, dihedrals, blocks, sidechains)
# Metadata(atoms=(...), ss=(...), residues=(...), dihedrals=(...), blocks=(...), sidechains=(...))

# julia> Metadata()
# Metadata(atoms=AtomMetadata[], ss=SecondaryStructureMetadata[], residues=Residue[], dihedrals = Diehdral[], blocks=[], sidechains=[])
# ```
# """
mutable struct Metadata
    atoms::Union{Nothing,Vector{AtomMetadata}}
    residues::Union{Nothing,Vector{Residue}}
    chains::Union{Nothing,Vector{Residue}}
    stree::Union{Nothing,SpanningTree}
end

Metadata(; atoms=nothing, residues=nothing, chains=nothing, stree=nothing) = begin
    if atoms != nothing
        # connectivity table
        #   atom i is connected to j is j>i
        table = [filter(x->x > at.index, at.connects) for at in atoms]
        stree = SpanningTree(table)
    end
    Metadata(atoms, residues, chains, stree)
end

Base.show(io::IO, b::Metadata) = begin
    print(io, "Metadata")
    for p in (:atoms, :residues, :chains)
        v = getproperty(b,p)
        s = v == nothing ? "n.d." : " list with $(length(v)) items"
        print(io, "\n   $(String(p)) = $s")
    end
end

# Base.@kwdef mutable struct Metadata

#     mst::Union{Nothing,Vector{Vector{Int}}} = nothing
#     atoms::Vector{AtomMetadata}            = Vector{AtomMetadata}()
#     residues::Vector{Residue}              = Vector{Residue}()
#     # chains::Vector{Chain}

#     #ss::Vector{SecondaryStructureMetadata} = Vector{SecondaryStructureMetadata}()
#     #dihedrals::Vector{Dihedral}            = Vector{Dihedral}()
#     #blocks::Vector{BlockMetadata}          = Vector{BlockMetadata}()
# end


# Base.show(io::IO, b::Metadata) = print(io, "Metadata(atoms=$(b.atoms), ss=$(b.ss), residues=$(b.residues), dihedrals=$(b.dihedrals))")


# @doc raw"""
#     renumber_residues!(atoms::Vector{AtomMetadata}[, start::Int64 = 1])

# Iterate over an array of AtomMetadata objects and renumber the list of `res_num` (starting at `start`).

# # Examples
# ```julia-repl
# julia> renumber_residues!(state.metadata.atoms)
# ```
# See also: [`AtomMetadata`](@ref)
# """
# function renumber_residues!(atoms::Vector{AtomMetadata}, start::Int64 = 1)

#     i_res_num::Int64 = atoms[1].res_num
#     for (index, atom) in enumerate(atoms)
#         atom.res_num = atom.res_num - i_res_num + start
#     end
# end