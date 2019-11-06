# -------------------------------------------------------------------------------------------------------------
#                                                DIHEDRAL

#Common.DTYPE.phi | Common.DTYPE.chi5
# module DIHEDRAL
#     @enum TYPE begin
#         phi   = -2
#         psi   = -1
#         omega =  0
#         chi1  =  1
#         chi2  =  2
#         chi3  =  3
#         chi4  =  4
#         chi5  =  5
#     end
# end


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
    #residue::Residue
    # type::DIHEDRAL.TYPE
    type::Enum{T<:Integer}
end
# Base.show(io::IO, b::Dihedral) = print(io, "Dihedral(a1=$(b.a1), a2=$(b.a2), a3=$(b.a3), a4=$(b.a4), movable=$(b.movable), residue=$(b.residue), type=$(b.dtype))")
function Base.show(io::IO, b::Dihedral)
    print(io, string(typeof(b)))
    for p in fieldnames(typeof(b))
        val = getproperty(b,p)
        print(io, "\n   $(String(p)) = $(val)")
    end
end


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

using Base.Cartesian

function rotate!(xyz::Array{Float64, 2}, dihedral::Dihedral, angle::Float64, tree::SpanningTree)
    a2 = dihedral.a2
    a3 = dihedral.a3
    
    lastidx = getindex(tree, a3)

    #pivot = @view xyz[a3,:]
    #axis = @view(xyz[a2,:]) - pivot

    axis = [0.0, 0.0, 0.0]
    # ox = xyz[a3,1]
    # oy = xyz[a3,2]
    # oz = xyz[a3,3]
    # axis[1] = xyz[a2,1] - ox
    # axis[2] = xyz[a2,2] - oy
    # axis[3] = xyz[a2,3] - oz
    @nexprs 3 k -> o_k = xyz[a3,k]
    @nexprs 3 k -> axis[k] = xyz[a2,k] - o_k
    R = Aux.rotation_matrix_from_axis_angle(axis, angle)

    @inbounds for i=1:lastidx
        j = tree.indices[i]
        
        @nexprs 3 k -> x_k = xyz[j,k]-o_k
        # x_1 = xyz[j,1] - o_1
        # x_2 = xyz[j,2] - o_2
        # x_3 = xyz[j,3] - o_3
        @nexprs 3 k -> xyz[j,k] = R[k,1]*x_1 + R[k,2]*x_2 + R[k,3]*x_3 + o_k
        # xyz[j,1] = rmat[1,1]*x_1 + rmat[1,2]*x_2 + rmat[1,3]*x_3 + o_1
        # xyz[j,2] = rmat[2,1]*x_1 + rmat[2,2]*x_2 + rmat[2,3]*x_3 + o_2
        # xyz[j,3] = rmat[3,1]*x_1 + rmat[3,2]*x_2 + rmat[3,3]*x_3 + o_3


    end

end


function rotate!(xyz::Array{Float64, 2}, dihedral::Dihedral, 
    angle::Float64, tree::SpanningTree, R::Array{Float64, 2})
    a2 = dihedral.a2
    a3 = dihedral.a3
    
    axis = [0.0, 0.0, 0.0]
    @nexprs 3 k -> o_k = xyz[a3,k]
    @nexprs 3 k -> axis[k] = xyz[a2,k] - o_k
    R = Aux.rotation_matrix_from_axis_angle(axis, angle, R)
    
    lastidx = getindex(tree, a3)
    @inbounds for i=1:lastidx
        j = tree.indices[i]
    #lastidx = size(dihedral.movable,1)
    #@inbounds for i=1:lastidx
    #    j = dihedral.movable[i]
        @nexprs 3 k -> x_k = xyz[j,k]-o_k
        @nexprs 3 k -> xyz[j,k] = R[k,1]*x_1 + R[k,2]*x_2 + R[k,3]*x_3 + o_k
    end
    xyz
end

