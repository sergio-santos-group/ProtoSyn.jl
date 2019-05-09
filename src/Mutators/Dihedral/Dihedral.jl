module Dihedral

using ..Common
using ..Aux
using ..Mutators
using ..Abstract

@doc raw"""
    DihedralMutator(dihedrals::Vector{Common.Dihedral}, angle_sampler::Function, p_mut::Float64, step_size::Float64)

Holds all necessary parameters for the correct simulation of dihedral movements.

# Arguments
- `dihedrals::Vector{Common.Dihedral}`: List of dihedrals avaliable to be rotated in dihedral movements.
- `angle_sampler::Function`: Function responsible for defining the new angle for the dihedral. Should return a `Float64`.
- `p_mut::Float64`: Probability of rotation of each dihedral.
- `step_size::Float64`: Scalar that defines the amount of change resulting for a dihedral movement.

# Examples
```julia-repl
julia> Mutators.Diehdral.DihedralMutator(dihedrals, randn, 0.05, 0.25)
DihedralMutator(dihedrals=68, p_pmut=0.05, angle_sampler=randn, step_size=0.25)
```
See also: [`run!`](@ref)
"""
@Base.kwdef mutable struct MutatorConfig{F <: Function} <: Abstract.MutatorConfig
    dihedrals::Vector{Common.Dihedral}
    angle_sampler::F
    p_mut::Float64 = 0.0
    step_size::Float64 = 0.0
end # struct
Base.show(io::IO, b::MutatorConfig) = print(io, "MutatorConfig(dihedrals=$(length(b.dihedrals)), angle_sampler=$(string(b.angle_sampler)), p_mut=$(b.p_mut), step_size=$(b.step_size))")


@doc raw"""
    run!(state::Common.State, mutator::DihedralMutator)

Iterate over a list of [`Common.Dihedral`](@ref) (`dihedrals`) and perform dihedral movements on the current
[`Common.State`](@ref). The probability of each dihedral undergo movements is defined in the
[`DihedralMutator`](@ref).`p_mut`. The new angle is obtained from [`DihedralMutator`](@ref).`angle_sampler`, who should
return a `Float64` in radians.
After movement, the [`Common.State`](@ref) is updated with the new conformation.
Returns the number of rotations performed.

# Examples
```julia-repl
julia> Mutators.Dihedral.run!(state, mutator)
```
See also: [`Common.rotate_dihedral!`](@ref Common)
"""
function apply!(state::Common.State, mutator::MutatorConfig)
    
    count::Int64 = 0
    for dihedral in mutator.dihedrals
        if rand() < mutator.p_mut
            Common.rotate_dihedral!(state.xyz, dihedral, mutator.angle_sampler())
            count += 1
        end
    end
    return count
end

end # module