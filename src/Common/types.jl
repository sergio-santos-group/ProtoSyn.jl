

const Optional = Union{Nothing, T} where T
const OptVector = Union{Nothing, Vector{T}} where T

abstract type AbstractMetadata end
abstract type AbstractAtomMetadata <: AbstractMetadata end

abstract type AbstractContainerMetadata <: AbstractMetadata end
abstract type AbstractResidueMetadata   <: AbstractContainerMetadata end
abstract type AbstractMoleculeMetadata  <: AbstractContainerMetadata end

#region AtomMetadata
Base.@kwdef mutable struct AtomMetadata <: AbstractAtomMetadata
    index::Int     = -1
    name::String   = "X"
    symbol::String = "X"
    parent::Optional{AbstractResidueMetadata} = nothing
end

@customshow AtomMetadata
#endregion



#region ResidueMetadata
Base.@kwdef mutable struct ResidueMetadata <: AbstractResidueMetadata
    index::Int   = -1
    name::String = "UNK"
    parent::Optional{AbstractMoleculeMetadata} = nothing
    children::Optional{Vector{AtomMetadata}} = nothing
end

@customshow ResidueMetadata
#endregion ResidueMetadata


const AdjacencyMetadata = IdDict{AtomMetadata, Vector{AtomMetadata}}


# #region ConnectivityMetadata
# struct ConnectivityMetadata
#     #size::Int
#     # table::Vector{Vector{AtomMetadata}}
#     table::IdDict{AtomMetadata, Vector{AtomMetadata}}
#     # ConnectivityMetadata(table::Vector{Vector{AtomMetadata}}, offset=0) = begin
#     #     table = map(enumerate(table)) do nr
#     #         n,r = nr
#     #         filter(x->x > n+offset, r)
#     #     end
#     #     new(size(table, 1), table)
#     # end
#     # ConnectivityMetadata(atoms::Vector{AtomMetadata}) = begin
#     #     table = map(atoms) do at
#     #         if at.connections == nothing
#     #             return Int[]
#     #         end
#     #         return filter(x->x > at.index, at.connections)
#     #     end
#     #     size = size(atoms, 1)
#     #     new(size, table)
#     # end
# end

# @customshow ConnectivityMetadata
# #endregion


#region MoleculeMetadata
Base.@kwdef mutable struct MoleculeMetadata <: AbstractMoleculeMetadata
    index::Int   = -1
    name::String = "MOL"
    children::OptVector{ResidueMetadata} = nothing
    adjacency::Optional{AdjacencyMetadata} = nothing
    
end

@customshow MoleculeMetadata
#endregion



#region DihedralMetadata
Base.@kwdef mutable struct BlockMetadata <: AbstractMetadata
    a1::Int
    a2::Int
    a3::Int
    a4::Int
    residues::OptVector{ResidueMetadata}
    atoms::OptVector{AtomMetadata}
    type::Optional{Enum{<:Integer}}
end
@customshow BlockMetadata
#endregion

# const Metadata = Dict{String, Vector{AbstractMetadata}}
const Metadata = Dict{Symbol, Vector{AbstractMetadata}}

Base.insert!(container::MoleculeMetadata, index::Integer, item::ResidueMetadata) = begin
    isa(container.children, Nothing) && (container.children = [])
    insert!(container.children, index, item)
    println("insert!(mol, $index, res); must recompile metadata")
    container
end

Base.push!(container::ResidueMetadata, item::AtomMetadata) = begin
    isa(container.children, Nothing) && (container.children = [])
    push!(container.children, item)
    container
end

Base.push!(container::MoleculeMetadata, item::ResidueMetadata) = begin
    isa(container.children, Nothing) && (container.children = [])
    push!(container.children, item)
    println("push!(mol, res); must recompile metadata")
    container
end

Base.delete!(container::MoleculeMetadata, item::ResidueMetadata) = begin
    i = findfirst(x->x===item, container.children)
    if !isa(i, Nothing)
        deleteat!(container.children, i)
        println("delete!(mol, res); must recompile metadata")
    end
    container
end

# Base.delete!(container::ResidueMetadata, item::AtomMetadata) = begin
#     i = findfirst(x->x===item, container.children)
#     if !isa(i, Nothing)
#         deleteat!(container.children, i)
#     end
#     container
# end

# REMARK POLYRULE
struct Rule
    from::Regex = r"VAL:C"
    to::Regex = r"GLY:N"
    # disallow intra-residue polymerization
    from::Regex = r"[A-Z]{3}:C"
    to::Regex = r"[A-Z]{3}:N"
    intra = false
    # allow intra-residue polymerization
    intra = true


    from::Regex = r"[A-Z]{3}:C(ola|adeus)"
    to::Regex = r"[A-Z]{3}:N"
    
    from::Regex = r"[A-Z]{3}:(AF)|(BF)"
    to::Regex   = r"[A-Z]{3}:(AT)|(BT)"

end




mutable struct Residue
    index::Int
    name::String
    parent::MoleculeMetadata
    children::Vector{AtomMetadata}
    next::Vector{Residue}
    prev::Vector{Residue}
    rule::Rule
end

mutable struct Molecule
    index::Int
    name::String
    children::Vector{Residue}
end













#region Energy
mutable struct Energy{T <: AbstractFloat}
    components::Dict{Symbol, T}
    total::T
end
Energy{T}() where {T <: AbstractFloat} = Energy{T}(Dict(), zero(T))

Energy(components::Dict{Symbol, T}) where {T <: AbstractFloat} = begin
    energy = Energy{T}(Dict(), zero(T))
    for key in keys(components)
        value = components[key]
        energy.components[key] = components[key]
        energy.total += value
    end
    energy
end

Base.setproperty!(e::Energy, component::Symbol, value::T) where {T <: AbstractFloat} = begin
    if component != :total
        components = getfield(e, :components)
        components[component] = value
    end
    setfield!(e, :total, value)
end

Base.getproperty(e::Energy, component::Symbol) = begin
    if component == :total
        return getfield(e, :total)
    elseif component == :components
        return getfield(e, :components)
    end
    getfield(e, :components)[component]
end
#endregion



#------------------------------------------------------------------------------
# const Optional = Union{Nothing, T} where T
# const AbstractContainerMetadata = Vector{AtomMetadata}
# const AtomMetadata = Int


mutable struct IndicesIterator
    singletons::Optional{Vector{AtomMetadata}}
    containers::Optional{Vector{AbstractContainerMetadata}}
    
    block::Int
    offset::Int
    bsize::Int
    size::Int
    current::Optional{Vector{AtomMetadata}}
end

IndicesIterator(s::Optional{Vector{AtomMetadata}}, c::Optional{Vector{<:AbstractContainerMetadata}}) = begin

    if isa(s, Nothing)
        block = 1
        bsize = size_s = 0
    else
        block = 0
        bsize = size_s = size(s, 1)
        
    end

    if !isa(c, Nothing) && size(c,1) > 0
        size_c = sum(map(x->length(x.children), c))
    else
        size_c = 0
    end

    if block==1 && size_c > 0
        bsize = size(c[1].children, 1)
        container = c[1].children
    else
        container = s
    end
    IndicesIterator(s, c, block, 0, bsize, size_s+size_c, container)
end

# struct BlockMetadata <: AbstractContainerMetadata
#     atoms::Optional{Vector{AtomMetadata}}
#     residues::Optional{Vector{AbstractContainerMetadata}}
# end

getiterator(r::ResidueMetadata)  = IndicesIterator(r.children, nothing)
getiterator(m::MoleculeMetadata) = IndicesIterator(nothing, m.children)
# getiterator(b::BlockMetadata) = IndicesIterator(b.atoms, b.residues)


Base.iterate(it::IndicesIterator, state::Int=1) = begin
    if state > it.size
        return nothing
    end
    
    # println("state=$state local=$(state-it.offset) block=$(it.block) offset=$(it.offset) bsize=$(it.bsize)")
    if state-it.offset > it.bsize
        it.block += 1
        it.offset += it.bsize
        it.current = it.containers[it.block].children
        it.bsize = size(it.current, 1)
    end
    (it.current[state-it.offset], state+1)
end



#----------------------------------------------------------------------------------------
### METHODS

function readpdb(filename::String, metadata=true)# where {T<:AbstractFloat}
    
    xyz = Vector{Vector{Float64}}()
    
    if metadata
        id2atom = Dict{Int, AtomMetadata}()
        bonds = Dict{Int, Vector{Int}}()

        atoms = Vector{AtomMetadata}()
        residues = Vector{ResidueMetadata}()
        molecules = Vector{MoleculeMetadata}()

        molecule = MoleculeMetadata(
            index=1,
            name="mol1",
            children=[]
        )

        push!(molecules, molecule)
        prev_rindex = -1

    end

    open(filename, "r") do fin
        for line in eachline(fin)

            if startswith(line, "ATOM")
                # save coordinates for this atom
                push!(xyz, [0.1*parse(Float64, line[n:n+7]) for n=31:8:47])

                # if metadata is requested
                if metadata

                    # extract residue information
                    rindex = parse(Int, line[23:26])
                    if prev_rindex != rindex
                        prev_rindex = rindex
                        rname = string(strip(line[18:20]))
                        residue = ResidueMetadata(
                            index=rindex,
                            name=rname,
                            children=[],
                            parent=molecules[end]
                        )
                        push!(residues, residue)
                        push!(molecules[end].children, residue)
                    end

                    # extract atom information
                    atom = AtomMetadata(
                        index  = parse(Int, strip(line[6:11])),
                        name   = string(strip(line[14:16])),
                        symbol = string(strip(line[77:78])),
                        parent = residues[end]
                    )
                    id2atom[atom.index] = atom
                    push!(atoms, atom)
                    push!(residues[end].children, atom)
                end

            elseif metadata && startswith(line, "CONECT")
                line = strip(line)
                indices = map(s->parse(Int, s), [line[n:n+4] for n=7:5:length(line)])
                append!(get!(bonds, indices[1], []), indices[2:end])
                
            elseif metadata && startswith(line, "TER")
                molecule = MoleculeMetadata(
                    index=molecules[end].index+1,
                    name=string("mol", molecules[end].index+1),
                    children=[]
                )
                push!(molecules, molecule)
            end

        end # end for line

    end # open


    if metadata
        
        for molecule in molecules
            bnds = IdDict{AtomMetadata, Vector{AtomMetadata}}()
            for at in getiterator(molecule)
                i = at.index
                if haskey(bonds, i)
                    bnds[id2atom[i]] = [id2atom[j] for j in bonds[i]]
                end
            end
        
            if length(bnds) > 0
                molecule.bonds = ConnectivityMetadata(bnds)
            end
        end
        
        return xyz, Metadata(
            :atoms => atoms,
            :residues => residues,
            :molecules => molecules
        )
    end
    reshape(vcat(xyz...),:,3), nothing
end


# struct MoleculeState
#     size::Int
#     coords::Matrix{Float64}
#     Forces::Matrix{Float64}
#     metadata::MoleculeMetadata
# end

# struct State
#     children::Vector{MoleculeState}
#     energy::Energy
#     size:Int
# end

# function Ebonds(state::State)

#     for s in state.children
#         for bond in bonds

#         end
#     end
    
# end




