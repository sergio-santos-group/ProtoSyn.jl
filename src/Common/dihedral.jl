# -------------------------------------------------------------------------------------------------------------
#                                                DIHEDRAL

#Common.DTYPE.phi | Common.DTYPE.chi5
module DIHEDRAL
    @enum TYPE begin
        phi   = -2
        psi   = -1
        omega =  0
        chi1  =  1
        chi2  =  2
        chi3  =  3
        chi4  =  4
        chi5  =  5
    end
end


@doc raw"""
    Dihedral(a1::Int64, a2::Int64, a3::Int64, a4::Int64, movable::Array{Int64, 1}, residue::Union{Common.Residue, Int64}, dtype::DIHEDRAL.TYPE)

Define a dihedral.

# Arguments
- `a1::Int64, a2::Int64, a3::Int64, a4::Int64`: *global* atom indices.
- `movable::Array{Int64, 1}`: List of *global* atom indices that will be moved during the dihedral movement in *this* residue.
- `residue::Union{Residue, Int64}`: [`Residue`](@ref) object that this dihedral belongs to (or Int64 identifying it).
- `dtype::DIHEDRAL.TYPE`: DIHEDRAL.TYPE

# Examples
```julia-repl
julia> Dihedral(2, 3, 5, 7, [5, 6], Common.Residue([1, 2, 3, 4, 5, 6], (...), "A", coil), phi)
Dihedral(a1=2, a2=3, a3=5, a4=7, movable=[5, 6], residue=Common.Residue(atoms=[1, 2, 3, 4, 5, 6], next=V, name=A), type=phi)
```
See also: [`Mutators.DihedralMutator`](@ref Mutators)
"""
mutable struct Dihedral
    a1::Int64
    a2::Int64
    a3::Int64
    a4::Int64
    movable::Vector{Int64}
    residue::Residue
    dtype::DIHEDRAL.TYPE
end
Base.show(io::IO, b::Dihedral) = print(io, "Dihedral(a1=$(b.a1), a2=$(b.a2), a3=$(b.a3), a4=$(b.a4), movable=$(b.movable), residue=$(b.residue), type=$(b.dtype))")


@doc raw"""
    rotate_dihedral!(xyz::Array{Float64, 2}, dihedral::Dihedral, angle::Float64)

Perform a dihedral movement, adding the provided `angle` (in radians). If the `dihedral.dtype` is "PHI"
or "PSI" the `dihedral.residue.next` is also rotated and this is propagated recursively until the end of
the molecule. 

# Examples
```julia-repl
julia> Mutators.Dihedral.rotate_dihedral!(state.xyz, dihedral, π/2)
```
See also: [`Dihedral`](@ref Mutators)
"""
function rotate_dihedral!(xyz::Array{Float64, 2}, dihedral::Dihedral, angle::Float64)
    rotate_dihedral!(xyz, dihedral.a2, dihedral.a3, angle, dihedral.dtype, dihedral.movable, dihedral.residue)
end


@doc raw"""
    rotate_dihedral!(xyz::Array{Float64, 2}, a2::Int64, a3::Int64, angle::Float64, dtype::DIHEDRAL.TYPE, movable::Vector{Int64}[, residue::Union{Residue, Nothing} = nothing])

Base dihedral movement function. Especifies all arguments used in dihedral rotation movement. 

# Examples
```julia-repl
julia> Mutators.Dihedral.rotate_dihedral!(xyz, dihedral.a2, dihedral.a3, π/2, dihedral.dtype, dihedral.movable, dihedral.residue)
```
See also: [`Aux.rotation_matrix_from_axis_angle`](@ref) [`Dihedral`](@ref Mutators) [`Crankshaft`](@ref Mutators)
"""
@inline function rotate_dihedral!(
    xyz::Array{Float64, 2},
    a2::Int64,
    a3::Int64,
    angle::Float64,
    dtype::DIHEDRAL.TYPE,
    movable::Vector{Int64},
    residue::Union{Residue, Nothing}=nothing)

    pivot = xyz[a2, :]'
    axis  = xyz[a3, :] - pivot'
    
    # Define the rotation matrix based on the rotation axis and angle
    rmat = Aux.rotation_matrix_from_axis_angle(axis, angle)

    #Rotate movable atoms pertaining to this dihedral
    xyz[movable, :] = (rmat * (xyz[movable, :] .- pivot)')' .+ pivot

    # Rotate all downstream residues
    if dtype <= DIHEDRAL.omega && residue != nothing
        idxs = Vector{Int64}()
        while residue.next != nothing
            residue = residue.next
            append!(idxs, residue.atoms)
        end
        xyz[idxs, :] = (rmat * (xyz[idxs, :] .- pivot)')' .+ pivot
    end
end


@doc raw"""
    rotate_dihedral_to!(xyz::Array{Float64, 2}, dihedral::Dihedral, target_angle::Float64)

Perform a dihedral movement, adding the necessary angle to meet the provided `target_angle` (in radians).

# Examples
```julia-repl
julia> Mutators.Dihedral.rotate_dihedral_to!(state.xyz, dihedral, π/2)
```
See also: [`apply_ss!`](@ref)
"""
function rotate_dihedral_to!(xyz::Array{Float64, 2}, dihedral::Dihedral, target_angle::Float64)

    #Calculate the current angle
    current_angle = Aux.calc_dih_angle(xyz[dihedral.a1, :], xyz[dihedral.a2, :], xyz[dihedral.a3, :], xyz[dihedral.a4, :])
    # Calculate necessary displacement to target angles
    displacement = target_angle - current_angle
    # Apply displacement and rotate
    rotate_dihedral!(xyz, dihedral, displacement)
end