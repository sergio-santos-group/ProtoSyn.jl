module Crankshaft

using ..Common
using ..Aux

@doc raw"""
    CrankshaftMutator(dihedrals::Vector{Common.Dihedral}, angle_sampler::Function, p_mut::Float64, step_size::Float64)

Holds all necessary parameters for the correct simulation of crankshaft movements.

# Arguments
- `dihedrals::Vector{Common.Dihedral}`: List of dihedrals avaliable to be rotated in crankshaft movements.
- `angle_sampler::Function`: Function responsible for defining the rotation angle. Should return a `Float64`.
- `p_mut::Float64`: Probability of rotation of each pair of alpha carbons.
- `step_size::Float64`: Scalar that defines the amount of change resulting for a crankshaft movement.

# Examples
```julia-repl
julia> Mutators.Crankshaft.CrankshaftMutator(dihedrals, randn, 0.05, 0.25)
CrankshaftMutator(dihedrals=68, angle_sampler=randn, p_pmut=0.05, step_size=0.25)

julia> Mutators.Crankshaft.CrankshaftMutator(dihedrals, randn)
CrankshaftMutator(dihedrals=68, angle_sampler=randn, p_pmut=0.0, step_size=0.0)
```
See also: [`run!`](@ref)
"""
mutable struct CrankshaftMutator
    dihedrals::Vector{Common.Dihedral}
    angle_sampler::Function
    p_mut::Float64
    step_size::Float64

    function CrankshaftMutator(dihedrals::Vector{Common.Dihedral}, angle_sampler::Function, p_mut::Float64, step_size::Float64)
        for dihedral in dihedrals
            dihedral.dtype == Common.DIHEDRAL.phi ? nothing : error("Tried to add a non-PHI dihedral to CrankshaftMutator ($dihedral)")
        end
        new(dihedrals, angle_sampler, p_mut, step_size)
    end
end
CrankshaftMutator(dihedrals::Vector{Common.Dihedral}, angle_sampler::Function; p_mut = 0.0, step_size = 0.0) = CrankshaftMutator(dihedrals, angle_sampler, p_mut, step_size)
Base.show(io::IO, b::CrankshaftMutator) = print(io, "CrankshaftMutator(dihedrals=$(length(b.dihedrals)), angle_sampler=$(string(b.angle_sampler)), p_mut=$(b.p_mut), step_size=$(b.step_size))")


@doc raw"""
    run!(state::Common.State, mutator::CrankshaftMutator)

Iterate over a list of [`Common.Dihedral`](@ref) and perform crankshaft movements on the current
[`Common.State`](@ref). The probability of each pair of alpha carbons undergo movement is defined in the
[`CrankshaftMutator`](@ref).`p_mut`. The new angle is obtained from [`CrankshaftMutator`](@ref).`angle_sampler`, who should
return a `Float64` in radians.
After movement, the [`Common.State`](@ref) is updated with the new conformation.
Returns the number of rotations performed.

# Examples
```julia-repl
julia> Mutators.Crankshaft.run!(state, mutator)
```
See also: [`rotate_crankshaft!`](@ref)
"""
@inline function run!(state::Common.State, mutator::CrankshaftMutator)

    l = length(mutator.dihedrals)
    count::Int64 = 0
    for i in 1:(l-1)
        for j in (i + 1):l
            if rand() < mutator.p_mut
                rotate_crankshaft!(state.xyz, mutator.dihedrals[i], mutator.dihedrals[j], mutator.angle_sampler())
                count += 1
            end
        end
    end
    return count
end


@doc raw"""
    rotate_crankshaft!(xyz::Array{Float64, 2}, dihedral1::Common.Dihedral, dihedral2::Common.Dihedral, angle::Float64)

Perform a crankshaft movement, adding the provided `angle` (in radians) to all atoms between the alpha carbon in `dihedral1` and `dihedral2`.
The last residue side_chain is also rotated and should be relaxed back to a equilibrium state (See: [`SteepestDescent`](@ref Drivers)).

# Examples
```julia-repl
julia> Mutators.Dihedral.rotate_crankshaft!(state.xyz, dihedral1, dihedral2, π/2)
```
"""
function rotate_crankshaft!(xyz::Array{Float64, 2}, dihedral1::Common.Dihedral, dihedral2::Common.Dihedral, angle::Float64)

    next = dihedral2.residue.next 
    dihedral2.residue.next = nothing
    Common.rotate_dihedral!(xyz, dihedral1.a3, dihedral2.a3, angle, dihedral1.dtype, dihedral1.movable, dihedral1.residue)
    dihedral2.residue.next = next
    
    Common.rotate_dihedral!(xyz, dihedral1.a3, dihedral2.a3, -angle, dihedral1.dtype, dihedral2.movable)
end

end