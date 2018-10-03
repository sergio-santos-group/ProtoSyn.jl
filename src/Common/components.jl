abstract type AbstractEnergy end

@doc raw"""
    NullEnergy()

Empty placeholder energy container.

# Examples
```julia-repl
julia> Common.NullEnergy()
Null
```
"""
struct NullEnergy <: AbstractEnergy end
Base.show(io::IO, b::NullEnergy) = print(io, "Null")

# ----------------------------------------------------------------------------------------------------------

@doc raw"""
    State(size::Int64, energy::AbstractEnergy, xyz::Array{Float64, 2}, forces::Array{Float64, 2}, atnames::Array{String, 1})

Define the current state of the system, containing the atoms positions, energy and forces applied.
If only `size::Int64` is provided, an empty State with the given `size` is created with zeros.

# Arguments
- `size::Int64`: Atom count in system.
- `energy::AbstractEnergy`: Current energy of the system (kJ mol⁻¹).
- `xyz::Array{Float64, 2}`: Atom positions in 3 dimensions.
- `forces::Array{Float64, 2}`: Forces applied in each dimension to each atom (kJ mol⁻¹ nm⁻¹)
- `atnames::Array{String, 1}`: List of atom names.

# Examples
```julia-repl
julia> Common.State(3)
Common.State(size=3, energy=Null, xyz=[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], forces=[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], atnames=String[])

julia> Common.State(2, Common.NullEnergy(), [1.1 1.1 1.1; 2.2 2.2 2.2], zeros(2, 3), ["C", "O"])
Common.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], atnames=["C", "O"])
```
"""
mutable struct State
    
    size::Int64
    energy::AbstractEnergy
    xyz::Array{Float64, 2}
    forces::Array{Float64, 2} # kJ mol⁻¹ nm⁻¹
    atnames::Array{String, 1}

end
State(n::Int64) = State(n, NullEnergy(), zeros(n, 3), zeros(n, 3), Array{String, 1}())
Base.show(io::IO, b::State) = print(io, "Common.State(size=$(b.size), energy=$(b.energy), xyz=$(b.xyz), forces=$(b.forces), atnames=$(b.atnames))")

# ----------------------------------------------------------------------------------------------------------

@doc raw"""
    Residue(atoms::Array{Int64, 1}, next::Union{Residue, Int64, Nothing}, name::String)

Define a residue as part of the system. 

# Arguments
- `atoms::Array{Int64, 1}`: list of atom *global* indices in this residue.
- `next::Union{Residue, Int64, Nothing}`: Next residue in the system, attached to this. Is preferably a Residue instance, but can in certain cases be the index in a Residue list or empty (Nothing).
- `name::String`: Name of the residue. If the residue describes an aminoacid, the correspondent letter is suggested.

# Examples
```julia-repl
julia> Common.Residue([1, 2, 3, 4], Common.Residue([5, 6, 7, 8], nothing, "V"), "E")
Common.Residue(atoms=[1, 2, 3, 4], next=V, name=E)
```
"""
mutable struct Residue
    atoms::Array{Int64, 1}
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

