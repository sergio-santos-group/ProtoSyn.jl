using Printf
using ProtoSyn
using ProtoSyn.Units

# Dict built when loading the appropriate LGrammar.
# Dihedral angles not in the dict can still be rotated, by manually defining
# the controling atom. For example, ALA-CB dihedral can be rotated by
# setting the dihedral value Δϕ on any of its children. By definition, 
# dihedral angles which only rotate the position of hydrogens are not
# considered in this dict.
chi_dict = Dict{String, Vector{String}}()

include("secondary-structure.jl")

abstract type DihedralType end

# --- Backbone dihedral angles

struct Phi <: DihedralType end; const phi = Phi()
struct Psi <: DihedralType end; const psi = Psi()
struct Omega <: DihedralType end; const omega = Omega()

(phi::Phi)(residue::Residue)     = begin
    if residue.parent == ProtoSyn.root(residue)
        @warn "Residue $residue has no phi angle"
        return nothing
    end
    return residue["C"]
end

(psi::Psi)(residue::Residue)     = begin
    if length(residue.children) == 0 
        @warn "Residue $residue has no psi angle"
        return nothing
    end
    return residue.children[1]["N"]
end

(omega::Omega)(residue::Residue) = begin
    if residue.parent == ProtoSyn.root(residue)
        @warn "Residue $residue has no omega angle"
        return nothing
    end
    return residue["CA"]
end