# ----------------------------------------------------------------------------------------------------------
#                                                RESIDUE

@doc raw"""
    Residue(atoms::Array{Int64, 1}, next::Union{Residue, Int64, Nothing}, name::String, ss::SecondaryStructureType)

Define a residue as part of the system.

# Arguments
- `atoms::Vector{Int64}`: list of atom *global* indices in this residue.
- `next::Union{Residue, Int64, Nothing}`: Next residue in the system, attached to this. Is preferably a Residue instance, but can in certain cases be the index in a Residue list or empty (Nothing).
- `name::String`: Name of the residue. If the residue describes an aminoacid, the correspondent letter is suggested.

# Examples
```julia-repl
julia> Common.Residue([1, 2, 3, 4], Common.Residue([5, 6, 7, 8], nothing, "V", Common.coil), "E")
Common.Residue(atoms=[1, 2, 3, 4], next=V, name=E)
```
"""
mutable struct Residue
    atoms::Vector{Int64}
    next::Union{Residue, Int64, Nothing}
    name::String
end
function Base.show(io::IO, b::Residue)
    if b.next == nothing
        print(io, "Common.Residue(atoms=$(b.atoms), next=nothing, name=$(b.name))")
    else
        print(io, "Common.Residue(atoms=$(b.atoms), next=$(b.next.name), name=$(b.name))")
    end
end