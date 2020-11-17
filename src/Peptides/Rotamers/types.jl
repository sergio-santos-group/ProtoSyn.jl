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
        chis *= " | $(titlecase(name)): ($(rad2deg(value))° ± $(rad2deg(sd)))"
    end
    print(io, "Rotamer{$T}: $(r.name)$chis")
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