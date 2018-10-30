# ----------------------------------------------------------------------------------------------------------
#                                                 STATE

@doc raw"""
    State(size::Int64, energy::AbstractEnergy, xyz::Array{Float64, 2}, forces::Array{Float64, 2}, metadata::Vector{Metadata})

Define the current state of the system, containing the atoms positions, energy and forces applied.
If only `size::Int64` is provided, an empty State with the given `size` is created with zeros.

# Arguments
- `size::Int64`: Atom count in system.
- `energy::AbstractEnergy`: Current energy of the system (kJ mol⁻¹).
- `xyz::Array{Float64, 2}`: Atom positions in 3 dimensions.
- `forces::Array{Float64, 2}`: Forces applied in each dimension to each atom (kJ mol⁻¹ nm⁻¹)
- `metadata::Metadata`: Metadata object.

# Examples
```julia-repl
julia> Common.State(3)
Common.State(size=3, energy=Null, xyz=[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], forces=[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], metadata=(...))

julia> Common.State(2, Common.NullEnergy(), [1.1 1.1 1.1; 2.2 2.2 2.2], zeros(2, 3), metadata)
Common.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], metadata=(...))
```
See also: [`Metadata`](@ref)
"""
mutable struct State

    size::Int64
    energy::AbstractEnergy
    xyz::Array{Float64, 2}
    forces::Array{Float64, 2} # kJ mol⁻¹ nm⁻¹

end
State(n::Int64) = State(n, NullEnergy(), zeros(n, 3), zeros(n, 3), Metadata())
Base.show(io::IO, b::State) = print(io, "State(size=$(b.size), energy=$(b.energy), xyz=$(b.xyz), forces=$(b.forces))")

function Base.iterate(st::State, idx = 1)

    if idx > size(st.xyz, 1)
        return nothing
    else
        return ((st.xyz[idx, :], st.metadata.atoms[idx]), idx+1)
    end
end