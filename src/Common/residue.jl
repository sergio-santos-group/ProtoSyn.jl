# ----------------------------------------------------------------------------------------------------------
#                                                RESIDUE

@doc raw"""
    SSTYPE

Enum: holds information regarding the secondary structure of each residue in the simulation.
"""
@enum SSTYPE begin
    coil  = 0
    alpha = 1
    beta  = 2
end


@doc raw"""
    Residue(atoms::Array{Int64, 1}, next::Union{Residue, Int64, Nothing}, name::String, ss::SecondaryStructureType)

Define a residue as part of the system.

# Arguments
- `atoms::Vector{Int64}`: list of atom *global* indices in this residue.
- `next::Union{Residue, Int64, Nothing}`: Next residue in the system, attached to this. Is preferably a Residue instance, but can in certain cases be the index in a Residue list or empty (Nothing).
- `name::String`: Name of the residue. If the residue describes an aminoacid, the correspondent letter is suggested.
- `ss::SSTYPE`: The intitial secondary structure type [`SSTYPE`](@ref) of this residue.

# Examples
```julia-repl
julia> Common.Residue([1, 2, 3, 4], Common.Residue([5, 6, 7, 8], nothing, "V", Common.coil), "E", Common.alpha)
Common.Residue(atoms=[1, 2, 3, 4], next=V, name=E, ss=alpha)
```
"""
mutable struct Residue
    atoms::Vector{Int64}
    next::Union{Residue, Int64, Nothing}
    name::String
    ss::SSTYPE
end
function Base.show(io::IO, b::Residue)
    if b.next == nothing
        print(io, "Common.Residue(atoms=$(b.atoms), next=nothing, name=$(b.name), ss=$(b.ss))")
    else
        print(io, "Common.Residue(atoms=$(b.atoms), next=$(b.next.name), name=$(b.name), ss=$(b.ss))")
    end
end