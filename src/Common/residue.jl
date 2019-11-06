# ----------------------------------------------------------------------------------------------------------
#                                                RESIDUE

#Common.SS.HELIX | Common.SS.SHEET
# module SS
# @enum TYPE begin
#     COIL  = 0
#     HELIX = 1
#     SHEET = 2
# end
# end


@doc raw"""
    Residue(atoms::Array{Int64, 1}, next::Union{Residue, Int64, Nothing}, name::String, ss::SecondaryStructureType)

Define a residue as part of the system.

# Arguments
- `atoms::Vector{Int64}`: list of atom *global* indices in this residue.
- `next::Union{Residue, Int64, Nothing}`: Next residue in the system, attached to this. Is preferably a Residue instance, but can in certain cases be the index in a Residue list or empty (Nothing).
- `name::String`: Name of the residue. If the residue describes an aminoacid, the correspondent letter is suggested.
- `ss::SS.TYPE`: Type of secondary structure

# Examples
```julia-repl
julia> Common.Residue([1, 2, 3, 4], Common.Residue([5, 6, 7, 8], nothing, "V", Common.coil), "E")
Common.Residue(atoms=[1, 2, 3, 4], next=V, name=E)
```
"""
Base.@kwdef mutable struct Residue
    #atoms::Vector{Int64}
    #next::Union{Residue, Int64, Nothing}
    index::Int = -1
    name::String = "UNK"
    ss::SS.TYPE = SS.COIL
end
function Base.show(io::IO, b::Residue)
    print(io, string(typeof(b)))
    for p in fieldnames(typeof(b))
        val = getproperty(b,p)
        print(io, "\n   $(String(p)) = $(val!=nothing ? val : "nothing")")
    end
end


# function Base.show(io::IO, b::Residue)
#     if b.next == nothing
#         print(io, "Common.Residue(atoms=$(b.atoms), next=nothing, name=$(b.name), ss=$(b.ss))")
#     else
#         print(io, "Common.Residue(atoms=$(b.atoms), next=$(b.next.name), name=$(b.name), ss=$(b.ss))")
#     end
# end