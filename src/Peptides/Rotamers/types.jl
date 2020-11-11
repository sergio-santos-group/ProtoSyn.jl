using ProtoSyn: State, setdihedral!, Residue
using ProtoSyn.Peptides.Dihedral: DihedralType
using Random

abstract type RotamerLibrary end

mutable struct Rotamer{T <: AbstractFloat}
    name::String
    probability::T
    chis::Dict{DihedralType, Tuple{T, T}}
end

function Base.show(io::IO, r::Rotamer{T}) where {T <: AbstractFloat}
    chis = ""
    for (index, (chi, (value, sd))) in enumerate(r.chis)
        chis *= "| Chi$index: ($value ± $sd)"
    end
    print(io, "Rotamer{$T}: $(r.name) - $(r.probability) $chis")
end

function apply!(state::State, rotamer::Rotamer, residue::Residue)
    @assert rotamer.name === residue.name "Tried to apply a $(rotamer.name) rotamer to residue $(residue.name)."

    for (chi, (value, sd)) in rotamer.chis
        _value = (randn() * sd) + value
        setdihedral!(state, chi(residue), _value)
    end
end

# ---

mutable struct BBD_RotamerLibrary{T <: AbstractFloat} <: RotamerLibrary
    name::String
    phis::Vector{T}
    psis::Vector{T}
    rotamers::Matrix{Vector{Rotamer}}
end

function Base.show(io::IO, r::BBD_RotamerLibrary)
    println(io, "Name: $(r.name) | Shape: $(size(r.rotamers))")
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
    return rl.rotamers[phi_index, psi_index]
end