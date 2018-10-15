#TODO: Update documentation
@doc raw"""
rotate_dihedral!(xyz::Array{Float64, 2}, dihedral::NewDihedral, angle::Float64)

Perform a dihedral movement, adding the provided `angle` (in radians). If the `dihedral.dtype` is "PHI"
or "PSI" the `dihedral.residue.next` is also rotated and this is propagated recursively until the end of
the molecule. 

# Examples
```julia-repl
julia> Mutators.Diehdral.rotate_dihedral(state.xyz, dihedral, π/2)
```
See also: [`run!`](@ref) [`Aux.rotation_matrix_from_axis_angle`](@ref)
"""
function rotate_dihedral!(
    xyz::Array{Float64, 2},
    dihedral::Dihedral,
    angle::Float64)

    rotate_dihedral!(xyz, dihedral.a2, dihedral.a3, angle, dihedral.dtype, dihedral.movable, dihedral.residue)
end

#TODO: Add documentation
@inline function rotate_dihedral!(
    xyz::Array{Float64, 2},
    a2::Int64,
    a3::Int64,
    angle::Float64,
    dtype::DIHEDRALTYPE,
    movable::Vector{Int64},
    residue::Union{Residue, Nothing}=nothing)

    pivot = xyz[a2, :]'
    axis  = xyz[a3, :] - pivot'
    
    # Define the rotation matrix based on the rotation axis and angle
    rmat = Aux.rotation_matrix_from_axis_angle(axis, angle)

    #Rotate movable atoms pertaining to this dihedral
    xyz[movable, :] = (rmat * (xyz[movable, :] .- pivot)')' .+ pivot

    
    # Rotate all downstream residues
    if (dtype < omega) && (residue != nothing)
        idxs = Vector{Int64}()
        while residue.next != nothing
            residue = residue.next
            append!(idxs, residue.atoms)
        end
        xyz[idxs, :] = (rmat * (xyz[idxs, :] .- pivot)')' .+ pivot
    end
end