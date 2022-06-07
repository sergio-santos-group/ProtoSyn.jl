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