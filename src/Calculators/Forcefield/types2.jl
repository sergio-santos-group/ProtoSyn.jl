using Base: @kwdef

export IForcefieldComponentType, Forcefield, genkey
export Atom, HarmonicBond, HarmonicAngle, CosineDihedral, NCosineDihedral
export AtomType, HarmonicBondType, HarmonicAngleType, CosineDihedralType, NCosineDihedralType


abstract type IForcefieldComponentType end
abstract type IForcefieldComponent end


@kwdef struct ForcefieldSpec
    name::String
    genpairs::Bool
    fudgeLJ::Float64
    fudgeQQ::Float64
    exclusion_depth::Int
    components::Dict{Symbol, Dict{String, IForcefieldComponentType}} = Dict()
end



#region ForcefieldTypes
mutable struct AtomType <: IForcefieldComponentType
    name::String
    z::Int
    m::Float64
    σ::Float64
    ϵ::Float64
    AtomType(n::String, z::Int, m::Float64, σ::Float64, ϵ::Float64) = begin
       new(n, z, m, 0.5*σ, sqrt(ϵ))
    end
end

mutable struct HarmonicBondType <: IForcefieldComponentType
    name::String
    r::Float64
    k::Float64
end

mutable struct HarmonicAngleType <: IForcefieldComponentType
    name::String
    θ::Float64
    k::Float64
end

mutable struct CosineDihedralType <: IForcefieldComponentType
    name::String
    θ::Float64
    k::Float64
    n::Float64
end

mutable struct NCosineDihedralType <: IForcefieldComponentType
    name::String
    θ::Vector{Float64}
    k::Vector{Float64}
    n::Vector{Float64}
end

#endregion ForcefieldTypes




#region ForcefieldComponents
mutable struct Atom <: IForcefieldComponent
    q::Float64
    type::AtomType
    exclusions::Opt{Vector{Int}}
    #Atom(q::Float64, t::AtomType, ex::Opt{Vector{Int}}) = begin
    #    new(sqrt(138.935458)*q, t, ex)
    #end
end
Atom(q::Float64, t::AtomType) = Atom(q, t, nothing)

mutable struct HarmonicBond <: IForcefieldComponent
    a1::Int
    a2::Int
    type::HarmonicBondType
end

mutable struct HarmonicAngle <: IForcefieldComponent
    a1::Int
    a2::Int
    a3::Int
    type::HarmonicAngleType
end

mutable struct CosineDihedral <: IForcefieldComponent
    a1::Int
    a2::Int
    a3::Int
    a4::Int
    type::CosineDihedralType
end

mutable struct NCosineDihedral <: IForcefieldComponent
    a1::Int
    a2::Int
    a3::Int
    a4::Int
    type::NCosineDihedralType
end

#endregion ForcefieldComponents






#region FFComponentContainers

abstract type IForcefieldConfig end

# mutable struct AtomContainer
#     atoms::Vector{Atom}
#     config::Opt{<:IForcefieldConfig}
# end
# AtomContainer(atoms::Vector{Atom}) = AtomContainer(atoms, nothing)

# AtomContainer() = AtomContainer(Atom[], nothing)

# Base.push!(c::AtomContainer, at::Atom) = begin
#     push!(c.atoms, at)
#     c
# end

mutable struct ComponentContainer{T <: IForcefieldComponent}
    items::Vector{T}
    config::Opt{<:IForcefieldConfig}
end

function ComponentContainer{T}() where {T <: IForcefieldComponent}
    ComponentContainer(T[], nothing)
end

Base.push!(c::ComponentContainer{T}, t::T) where {T <: IForcefieldComponent} = begin
    push!(c.items, t)
    c
end



export Topology

#const ComponentContainer = Union{AtomContainer, Vector{<:IForcefieldComponent}}

#const Topology = Dict{Symbol, ComponentContainer}
struct Topology
    size::Int
    components::Dict{Symbol, ComponentContainer}
end
Topology() = Topology(-1, Dict())
Topology(n::Int) = Topology(n, Dict())

Base.push!(top::Topology, item::IForcefieldComponent) = begin
    type = typeof(item)
    key = nameof(type)
    push!(
        # get!(top, key, key==:Atom ? AtomContainer() : Vector{type}() ),
        #get!(top.components, key, key==:Atom ? AtomContainer() : Vector{type}() ),
        get!(top.components, key, ComponentContainer{type}()),
        item
    )
    top
end



# +(n::Int, item::T) where {T<:IForcefieldComponent} = begin
#     indices = [
#         n+getfield(item, Symbol(:a,i))
#         for i = 1:fieldcount(T)-1]
#             T(indices..., item.type)
#         end
        
#         +(item::T, n::Int) where {T<:IForcefieldComponent} = n + item
        
import Base.+
+(n::Int, t::Atom) = begin
    at = Atom(t.q, t.type)
    if t.exclusions !== nothing
        at.exclusions = n .+ t.exclusions
        at.exclusions[end] = -1
    end
    at
end
+(n::Int, t::HarmonicBond) = HarmonicBond(t.a1+n, t.a2+n, t.type)
+(n::Int, t::HarmonicAngle) = HarmonicAngle(t.a1+n, t.a2+n, t.a3+n, t.type)
+(n::Int, t::CosineDihedral) = CosineDihedral(t.a1+n, t.a2+n, t.a3+n, t.a4+n, t.type)
+(n::Int, t::NCosineDihedral) = NCosineDihedral(t.a1+n, t.a2+n, t.a3+n, t.a4+n, t.type)

+(t::T, n::Int) where {T <: IForcefieldComponent} = n+t

# +(t::Atom, n::Int) = n+t
# +(t::HarmonicBond, n::Int) = n+t
# +(t::HarmonicAngle, n::Int) = n+t
# +(t::CosineDihedral, n::Int) = n+t
# +(t::NCosineDihedral, n::Int) = n+t


import Base.*
*(n::Int, top::Topology) = begin
    tp = Topology(n * top.size)
    for i in 0:n-1
        offset = i * top.size
        for container in values(top.components)
            foreach(c->push!(tp,offset+c), container.items)
        end
    end
    tp
end
*(top::Topology, n::Int) = n * top

+(t1::Topology, t2::Topology) = begin
    t12 = Topology(t1.size + t2.size)
    for container in values(t1.components)
        foreach(c->push!(t12, 0+c), container.items)
    end
    offset = t1.size
    for container in values(t2.components)
        foreach(c->push!(t12, offset+c), container.items)
    end
    t12
end


using Printf: @sprintf

function Base.write(io::IO, top::Topology, st::State)
    if !haskey(top.components, :Atom)
        return write(io, st)
    end
    
    x = st.coords
    atoms = top.components[:Atom].items
    println(io, "MODEL")
    for i = 1:st.size
        println(io, @sprintf("ATOM %6d %4s UNK A%4d    %8.3f%8.3f%8.3f",
            i, atoms[i].type.name, i, 10.0*x[i,1], 10.0*x[i,2], 10.0*x[i,3])
        )
    end
    println(io, "ENDMDL")
end


#     c = fieldcount(T)
#     indices = [n+getfield(item, Symbol(:a,i)) for i = 1:c-1]
#     T(indices..., item.type)
# end

# +(t1::Topology, t2::Topology) = begin
#     t = Topology()
#     for
# end
#endregion FFComponentContainers



function aggregate(components::ComponentContainer...)
    function(state::State, forces::Bool)
        forces && fill!(state.forces, 0.0)
        e = 0.0
        for comp in components
            e += Calculators.eval!(state, comp, forces)
        end
        state.energy.total = e
    end
end

aggregate(top::Topology) = aggregate(values(top.components)...)

