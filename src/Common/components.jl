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

@enum SecondaryStructureType begin
    coil  = 0
    alpha = 1
    beta  = 2
end


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
    ss::SecondaryStructureType
end
function Base.show(io::IO, b::Residue)
    if b.next == nothing
        print(io, "Common.Residue(atoms=$(b.atoms), next=nothing, name=$(b.name))")
    else
        print(io, "Common.Residue(atoms=$(b.atoms), next=$(b.next.name), name=$(b.name))")
    end
    
end

# -------------------------------------------------------------------------------------------------------------

@enum DIHEDRALTYPE begin
    phi   = 0
    psi   = 1
    omega = 2
    chi1  = 3
    chi2  = 4
    chi3  = 5
    chi4  = 6
    chi5  = 7
end

@doc raw"""
    MutableDihedral(a1::Int64, a2::Int64, a3::Int64, a4::Int64, movable::Array{Int64, 1}, residue::Union{Common.Residue, Int64}, dtype::DIHEDRALTYPE)

Define a dihedral.

# Arguments
- `a1::Int64, a2::Int64, a3::Int64, a4::Int64`: *global* atom indices.
- `movable::Array{Int64, 1}`: List of *global* atom indices that will be moved during the dihedral movement in *this* residue.
- `residue::Union{Common.Residue, Int64}`: Residue that this dihedral belongs to. Should be a [`Common.Residue`](@ref) object.
- `dtype::String`: Dihedral type (i.e. "PHI", "PSI", ...)

# Examples
```julia-repl
julia> Mutators.Diehdral.NewDihedral(2, 3, 5, 7, [5, 6], Common.Residue([1, 2, 3, 4, 5, 6], (...), "A"), "PSI")
Dihedral(a1=2, a2=3, a3=5, a4=7, movable=[5, 6], residue=Common.Residue(atoms=[1, 2, 3, 4, 5, 6], next=V, name=A), type="PSI")
```
See also: [`load_topology`](@ref)
"""
mutable struct Dihedral
    a1::Int64
    a2::Int64
    a3::Int64
    a4::Int64
    movable::Vector{Int64}
    residue::Residue
    dtype::DIHEDRALTYPE
end
Base.show(io::IO, b::Dihedral) = print(io, "Common.Dihedral(a1=$(b.a1), a2=$(b.a2), a3=$(b.a3), a4=$(b.a4), movable=$(b.movable), residue=$(b.residue), type=$(b.dtype))")

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
    atnames::Vector{String}

end
State(n::Int64) = State(n, NullEnergy(), zeros(n, 3), zeros(n, 3), Array{String, 1}(), [])
Base.show(io::IO, b::State) = print(io, "Common.State(size=$(b.size), energy=$(b.energy), xyz=$(b.xyz), forces=$(b.forces), atnames=$(b.atnames))")

# ----------------------------------------------------------------------------------------------------------

#TODO: Document structure
mutable struct CallbackObject
    freq::Int64
    callback::Function
end