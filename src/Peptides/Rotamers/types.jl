using ProtoSyn: State, setdihedral!, Residue
using ProtoSyn.Peptides.Dihedral: DihedralType
using Random
using StatsBase

abstract type RotamerLibrary end

mutable struct Rotamer{T <: AbstractFloat}
    name::String
    chis::Dict{DihedralType, Tuple{T, T}}
end

function Base.show(io::IO, r::Rotamer{T}) where {T <: AbstractFloat}
    chis = ""
    for index in 1:length(r.chis)
        name = "chi$index"
        value, sd = r.chis[getfield(ProtoSyn.Peptides.Dihedral, Symbol(name))]
        chis *= " | $(titlecase(name)): ($(value)° ± $sd)"
    end
    print(io, "Rotamer{$T}: $(r.name)$chis")
end

function apply!(state::State, rotamer::Rotamer, residue::Residue)
    @assert rotamer.name === residue.name "Tried to apply a $(rotamer.name) rotamer to residue $(residue.name)."

    for (chi, (value, sd)) in rotamer.chis
        _value = (randn() * sd) + value
        setdihedral!(state, chi(residue), _value)
    end
end

# ---

mutable struct RotamerStack{T <: AbstractFloat}
    rotamers::Vector{Rotamer{T}}
    weights::Weights{T, T, Array{T, 1}}
end

RotamerStack(::Type{T}) where {T <: AbstractFloat} = begin
    RotamerStack(Vector{Rotamer{T}}(), Weights(Vector{T}()))
end

RotamerStack() = RotamerStack(ProtoSyn.Units.defaultFloat)

function sample(rs::RotamerStack{T}, n::Int = -1) where {T <: AbstractFloat}
    if n <= 0
        return StatsBase.sample(rs.rotamers, rs.weights)
    else
        return StatsBase.sample(rs.rotamers[1:n], rs.weights[1:n])
    end
end

Base.push!(rs::RotamerStack{T}, rotamer::Rotamer{T}, weight::Float64) where {T <: AbstractFloat} = begin

    push!(rs.rotamers, rotamer)
    rs.weights = Weights(vcat(rs.weights, Weights([weight])))
end

# ---

mutable struct BBD_RotamerLibrary{T <: AbstractFloat} <: RotamerLibrary
    name::String
    phis::Vector{T}
    psis::Vector{T}
    rotamer_stacks::Matrix{RotamerStack}
end

function Base.show(io::IO, r::BBD_RotamerLibrary)
    println(io, "Name: $(r.name) | Shape: $(size(r.rotamer_stacks))")
end

function findnearest(a, x)
    length(a) > 0 || return nothing
    r = searchsorted(a,x)
    length(r) > 0 && return first(r)
    last(r) < 1 && return first(searchsorted(a,a[first(r)]))
    first(r) > length(a) && return first(searchsorted(a,a[last(r)]))
    x-a[last(r)] <= a[first(r)]-x && return first(searchsorted(a,a[last(r)]))
    x-a[last(r)] > a[first(r)]-x && return first(searchsorted(a,a[first(r)]))
end


Base.getindex(rl::R, phi::T, psi::T) where {R <: RotamerLibrary, T <: AbstractFloat} = begin
    phi_index = findnearest(rl.phis, phi)
    psi_index = findnearest(rl.psis, psi)
    return rl.rotamer_stacks[phi_index, psi_index]
end