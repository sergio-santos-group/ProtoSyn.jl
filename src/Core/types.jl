using Base: @kwdef

export Atom, Residue, Molecule
export Link, LinkedResidue
export State, Dihedral, AxisRotatableBlock


# optional type
export Opt
const Opt = Union{Nothing, T} where T


# connectivity graph
export ConnectGraph
const ConnectGraph = Dict{Int, Vector{Int}}
const ConnectGraphByName = Dict{String, Vector{String}}


# ===================================================================

# @kwdef mutable struct _Atom{T}
#     id::Int
#     name::String
#     symbol::String
#     x::Float64
#     y::Float64
#     z::Float64
#     parent::Opt{T} = nothing
# end
mutable struct _Atom{T}
    id::Int
    symbol::String
    name::String
    x::Float64
    y::Float64
    z::Float64
    parent::Opt{T}
    
    _Atom{T}(id::Int, symbol::String, name::String, x::Float64, y::Float64, z::Float64) where {T}= begin
        new{T}(id, symbol, name, x, y, z, nothing)
    end

end

# ===================================================================

#@kwdef
mutable struct Residue
    name::String
    atoms::Vector{_Atom{Residue}}
    bonds::Opt{ConnectGraph}
    atomsbyname::Dict{String, _Atom{Residue}}
    Residue(name::String) = begin
        r = new()
        r.name = name
        r.atoms = []
        r.bonds = nothing
        r.atomsbyname = Dict{String, _Atom{Residue}}()
        r
    end
end
const Atom = _Atom{Residue}



Residue(name::String, atoms::Vector{Atom}) = begin
    r = Residue(name)
    push!(r, atoms...)
end


Residue(name::String, atoms::Vector{Atom}, bonds::ConnectGraph) = begin
    r = Residue(name)
    push!(r, atoms...)
    
    natoms = length(r.atoms)

    # add bonds. One could simply take the "bonds" object as is
    # and assign it to the residue. However, this way one ensures
    # that the connection graph is consistent and bi-directional.
    r.bonds = ConnectGraph()
    for (i,jj) in bonds
        @assert i <= natoms && r.atoms[i].id == i "atom with id=$(i) must exist"
        for j in jj
            @assert j <= natoms && r.atoms[j].id == j "atom with id=$(j) must exist"

            ai = i
            aj = j
            if ai>aj
                ai,aj=aj,ai
            end
            c = get!(r.bonds, ai, [])
            if aj ∉ c
                push!(c, aj)
                push!(get!(r.bonds, aj, []), ai)
            end
        end
    end
    r
end

Residue(name::String, atoms::Vector{Atom}, bonds::ConnectGraphByName) = begin
    name2atom = Dict{String, Atom}()
    for (index,atom) in enumerate(atoms)
        name2atom[atom.name] = atom
        atom.id = index
    end

    ibonds = ConnectGraph(
        name2atom[k].id => map(s -> name2atom[s].id, bonds[k])
        for k in keys(bonds)
    )
    Residue(name, atoms, ibonds)
end


# ===================================================================


export ResidueLib
const ResidueLib = Dict{String, Residue}


# ===================================================================

@kwdef mutable struct _Link{T}
    residue1::T
    residue2::T
    bonds::Opt{ConnectGraph} = nothing
end


# ===================================================================

@kwdef mutable struct LinkedResidue
    id::Int = -1
    source::Residue
    offset::Int = 0
    links::Vector{_Link{LinkedResidue}} = []
    flag::Bool = false
end
const Link = _Link{LinkedResidue}


# ===================================================================

@kwdef mutable struct Molecule
    name::String = "molecule"
    residues::Vector{LinkedResidue} = []
    links::Vector{Link} = []
    bonds::Opt{ConnectGraph} = nothing
    coherent::Bool = false
    offset::Int = 0
    size::Int = 0
end


# ===================================================================

#region Energy
abstract type  AbstractEnergy{T <: AbstractFloat} end

mutable struct Energy{T} <: AbstractEnergy{T}
    total::T
    components::Dict{Symbol, T}
end
Energy(e::T) where {T <: AbstractFloat} = Energy(e, Dict{Symbol, T}())
Energy{T}() where {T <: AbstractFloat} = Energy(T(NaN))

# Energy(components::Dict{Symbol, T}) where {T <: AbstractFloat} = begin
#     # energy = Energy(T(0))
    
#     # total = T(0)
#     # d = energy.components
    
#     # for key in keys(components)
#     #     value = components[key]
#     #     d[key] = components[key]
#     #     total += value
#     # end
    
#     # energy.total = total
#     # energy
#     Energy{T}(sum(values(components)), copy(components))
# end

Base.setproperty!(e::Energy{T}, k::Symbol, value::T) where {T <: AbstractFloat} = begin
    if k != :total
        components = getfield(e, :components)
        components[k] = value
    end
    setfield!(e, :total, value)
end

Base.getproperty(e::Energy{T}, k::Symbol) where {T <: AbstractFloat} = begin
    if k == :total
        return getfield(e, :total)
    elseif k == :components
        return getfield(e, :components)
    end
    getfield(e, :components)[k]
end

Base.copy(e::Energy) = Energy(e.total, copy(e.components))

Base.copy!(dst::Energy{T}, src::Energy{T}) where {T <: AbstractFloat} = begin
    copy!(dst.components, src.components)
    dst.total = src.total
    dst
end

# Base.convert(::Type{T}, x::Energy) where {T<:AbstractFloat} = Energy{T}(T(x.total), Dict{Symbol, T}(k=>T(v) for (k,v) in x.components))

#endregion

#--------------------------------------------------------------------

# mutable struct VerletList
#     capacity::Int               # max size of list 
#     cutoff::Float64             # cutoff (nm)
#     buffer::Float64             # buffer (nm)
#     offset::Opt{Vector{Int}}    # pointer particle->list position
#     list::Opt{Vector{Int}}      # nblist
#     # _
# end
# VerletList(n::Int) = VerletList(n, -1.0, 0.0, zeros(Int, n), zeros(Int, n))


# # function Base.copy!(dst::VerletList, src::VerletList)
# #     copy!(dst.pointer, src.pointer)
# #     copy!(dst.list, src.list)
# #     dst.cutoff = src.cutoff
# #     dst.buffer = src.buffer
# #     return dst
# # end




# Base.resize!(vlist::VerletList, n::Int) = begin
#     resize!(vlist.list, n)
#     vlist.capacity = n
#     vlist
# end

# function update!(vlist::VerletList, coords::Matrix{Float64})
#     if vlist.cutoff < 0.0
#         return vlist
#     end
    
#     ptr = 1
#     natoms = size(coords, 1)
#     cutsq = (vlist.cutoff + vlist.buffer)^2

#     for i = 1:natoms
#         vlist.offset[i] = ptr
#         for j = (i+1):natoms
#             @nexprs 3 u -> vij_u = coords[i,u] - coords[j,u]
#             if @dot(u, vij_u, vij_u) < cutsq
#                 vlist.list[ptr] = j
#                 ptr += 1
#                 if ptr == vlist.capacity
#                     resize!(vlist, vlist.capacity + 1000)
#                 end
#             end
#         end
#         vlist.list[ptr] = -1
#         ptr += 1

#         if (i < natoms) && (ptr == vlist.capacity)
#             resize!(vlist, vlist.capacity + 1000)
#         end

#     end
#     vlist
# end

# ===================================================================

mutable struct State{T <: AbstractFloat}
    size::Int
    coords::Array{T, 2}
    forces::Opt{Array{T, 2}}
    velocs::Opt{Array{T, 2}}
    energy::Energy{T}
    # pairlist::Opt{VerletList}
    # function State(s::State)
    #     state = new()
    #     state.size = s.size
    #     for field in fieldnames(State)[2:end]
    #         src = getproperty(s, field)
    #         if src !== nothing
    #             setproperty!(state, field, copy(src))
    #         end 
    #     end
    #     state
    # end
end

State{T}(size::Int) where {T<:AbstractFloat} = begin
    State(size, zeros(T, size, 3), nothing, nothing, Energy{T}())
end

# State(size::Int) = State(size, zeros(size,3), nothing, nothing, nothing, nothing)

function State(s::State{T}) where {T<:AbstractFloat}
    state = State{T}(s.size)
    for field in fieldnames(State)
    #for i = 2:fieldcount(State)
    #    field = fieldname(State,i)
        src = getproperty(s, field)
        if src !== nothing
            setproperty!(state, field, copy(src))
        end 
    end
    state
end

# export nbupdate!
# nbupdate!(s::State) = update!(s.pairlist, s.coords)



# macro gencopy!(t::Symbol, dst::Symbol, src::Symbol)
#     ex = quote end
#     type = getfield(@__MODULE__, t)
#     for comp in fieldnames(type)
#         push!(ex.args, quote
#             if $(dst).$(comp) !== nothing
#                 copy!($(dst).$(comp), $(src).$(comp))
#             end
#         end)
        
#     end
#     esc(ex)
# end



Base.copy!(x::Nothing, y::Nothing)::Nothing = nothing

Base.copy!(dst::State{T}, src::State{T}) where {T <: AbstractFloat} = begin
    @assert dst.size == src.size "copy! source and destiny state sizes must be equal"
    
    for i = 2:fieldcount(State)
        field = fieldname(State,i)
        copy!(getproperty(dst, field), getproperty(src, field))
    end

    # if (src.coords !== nothing) copy!(dst.coords, src.coords) end
    # if (src.forces !== nothing) copy!(dst.forces, src.forces) end
    # if (src.velocs !== nothing) copy!(dst.velocs, src.velocs) end
    # if (src.energy !== nothing) copy!(dst.energy, src.energy) end
    # if (src.pairlist !== nothing) copy!(dst.pairlist, src.pairlist) end

    dst
end






# ===================================================================

abstract type AbstractAxisRotatableBlock end

#@kwdef 
mutable struct Dihedral <: AbstractAxisRotatableBlock
    a0::Int
    a1::Int
    a2::Int
    a3::Int
    start::Int
    type::Opt{Enum{<:Integer}}
end
Dihedral(a0::Int, a1::Int, a2::Int, a3::Int) = Dihedral(a0, a1, a2, a3, a2, nothing)
Dihedral(a0::Int, a1::Int, a2::Int, a3::Int, type::Opt{Enum{<:Integer}}) = 
    Dihedral(a0, a1, a2, a3, a2, type)

#@kwdef
mutable struct AxisRotatableBlock <: AbstractAxisRotatableBlock
    a1::Int
    a2::Int
    start::Int
    #enpoints::Vector{Int} = []
    # movable::Opt{Vector{Int}} = nothing
    # terminals::Opt{Vector{Int}} = nothing
end



# ===================================================================




iterbyatom(r::Residue) = (a for a in r.atoms)
iterbyatom(lr::LinkedResidue) = (a for a in lr.source.atoms)
iterbyatom(m::Molecule) = (a for lr in m.residues for a in lr.source.atoms)

iterbyresidue(m::Molecule) = (lr for lr in m.residues)




# ===================================================================


function update!(mol::Molecule)
    
    # do nothing if this molecule is marked as coherent!!
    if mol.coherent
        return mol
    end

    # reset connectivity map
    if mol.bonds === nothing
       mol.bonds = ConnectGraph()
    else
        empty!(mol.bonds)
    end

    offset = 0

    # copy local (residue) bonds to molecule. This is a mere copy of the
    # ConnectGraph of the source with the due offset for this residue. Also,
    # update the offset of each l-residue
    for (rid, lr) in enumerate(mol.residues)
        for id1 in keys(lr.source.bonds)
            mol.bonds[id1 + offset] = offset .+ lr.source.bonds[id1]
        end
        lr.id = rid
        lr.offset = offset
        offset += length(lr.source)
    end
    
    # by now, the offset corresponds to the total number of
    # atoms in the molecule. Hence, one can set the size of the
    # molecule the same as the offset!
    mol.size = offset

    # add link bonds to the molecular ConnectGraph. One has to take into 
    # consideration the offsets of both l-residues.
    for link in mol.links
        offset1 = link.residue1.offset
        offset2 = link.residue2.offset
        
        for id1 in keys(link.bonds)
            # ConnectGraph keys (id1) are always local to residue1
            # whereas the values (id2) are local to residue2. Hence,
            # id1 and id2 need to be offset by offset1 and offset2, respectively
            bonds = get!(mol.bonds, id1, [])
            for id2 in link.bonds[id1]
                push!(get!(mol.bonds, id1+offset1, []), id2+offset2)
                push!(get!(mol.bonds, id2+offset2, []), id1+offset1)
            end
        end

    end

    # make sure to mark this molecule as being coherent
    mol.coherent = true

    # return the molecule
    mol
end


# ===================================================================


# function bind(f1::LinkedResidue, f2::LinkedResidue, rule::T) where {T<:Function}
function bind(rule::T, lr1::LinkedResidue, lr2::LinkedResidue, mol::Molecule) where {T<:Function}
    
    if (lr1 ∉ mol.residues)
        error("lr1 must be a child of the given molecule")
    elseif (lr2 ∉ mol.residues)
        error("lr2 must be a child of the given molecule")
    end

    link = rule(lr1, lr2)
    if link===nothing
        return
    end
    push!(lr1, link)
    push!(lr2, link)
    push!(mol, link)
    link
end