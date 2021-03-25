using ProtoSyn: State, setdihedral!, Residue
using ProtoSyn.Peptides.Dihedral: DihedralType
using Random
using StatsBase
using Printf

abstract type RotamerLibrary end

mutable struct Rotamer{T <: AbstractFloat}
    name::String
    chis::Dict{DihedralType, Tuple{T, T}}
end

function Base.show(io::IO, r::Rotamer{T}) where {T <: AbstractFloat}
    print(io, "Rotamer{$T}: $(as_string(r))\n")
end

function as_string(r::Rotamer{T}) where {T <: AbstractFloat}
    chis = ""
    for index in 1:length(r.chis)
        name = "chi$index"
        value, sd = r.chis[getfield(ProtoSyn.Peptides.Dihedral, Symbol(name))]
        chis *= @sprintf " | %s: %7.1f° ± %4.1f" titlecase(name) rad2deg(value) rad2deg(sd)
    end
    return "$(r.name)$chis"
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

function Base.show(io::IO, r::RotamerStack)
    println(io, "\n+"*repeat("-", 104)*"+")
    @printf(io, "| %-5s | %10s | %-80s |\n", "Index", "Probability", "Rotamer description")
    println(io, "+"*repeat("-", 104)*"+")
    for (index, rotamer) in enumerate(r.rotamers)
        @printf(io, "| %-5d | %10.2f%% | %-80s |\n", index, r.weights[index]*100, as_string(rotamer))
    end
    println(io, "+"*repeat("-", 104)*"+")
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

Base.copy(rot_lib::Rotamers.BBD_RotamerLibrary{T}) where {T <: AbstractFloat} = begin
    return BBD_RotamerLibrary(
        rot_lib.name,
        copy(rot_lib.phis),
        copy(rot_lib.psis),
        copy(rot_lib.rotamer_stacks)
    )
end