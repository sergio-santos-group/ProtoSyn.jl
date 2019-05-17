module Dihedral

using ..Common
using ..Aux
using ..Mutators
using ..Abstract

@doc raw"""
    MutatorConfig(; dihedrals::Vector{Common.Dihedral} = [], angle_sampler::Function = rand, p_mut::Float64 = 0.0, step_size::Float64 = 0.0) <: Abstract.MutatorConfig

Holds all necessary parameters for the correct simulation of dihedral movements.

# Arguments
- `dihedrals::Vector{Common.Dihedral}`: List of dihedrals avaliable to be rotated in dihedral movements (Default: empty).
- `angle_sampler::Function`: Function responsible for defining the new angle for the dihedral. Should return a `Float64` (Default: rand()).
- `p_mut::Float64`: Probability of rotation of each dihedral (Default: 0.0).
- `step_size::Float64`: Scalar that defines the amount of change resulting for a dihedral movement (Default: 0.0).

# Examples
```julia-repl
julia> Mutators.Dihedral.MutatorConfig(dihedrals = dihedrals, angle_sampler = randn, pmut = 0.05, step_size = 0.25)
Dihedral.MutatorConfig(dihedrals=68, angle_sampler=randn, pmut=0.05 step_size=0.25)

julia> Mutators.Dihedral.MutatorConfig()
Dihedral.MutatorConfig(dihedrals=0, angle_sampler=rand, p_pmut=0.0, step_size=0.0)
```
See also: [`apply!`](@ref)
"""
Base.@kwdef mutable struct MutatorConfig{F <: Function} <: Abstract.MutatorConfig
    
    # Paramters                        # Defaults
    dihedrals::Vector{Common.Dihedral} = Vector{Common.Dihedral}()
    angle_sampler::F                   = rand
    p_mut::Float64                     = 0.0
    step_size::Float64                 = 0.0
end # end struct
Base.show(io::IO, b::MutatorConfig) = print(io, "Dihedral.MutatorConfig(dihedrals=$(length(b.dihedrals)), angle_sampler=$(string(b.angle_sampler)), p_mut=$(b.p_mut), step_size=$(b.step_size))")


@doc raw"""
    run!(state::Common.State, mutator::MutatorConfig)

Iterate over a list of [`Common.Dihedral`](@ref) (`dihedrals`) and perform dihedral movements on the current
[`Common.State`](@ref). The probability of each dihedral undergo movements is defined in the
[`MutatorConfig`](@ref).`p_mut`. The new angle is obtained from [`MutatorConfig`](@ref).`angle_sampler`, who should
return a `Float64` in radians.
After movement, the [`Common.State`](@ref) is updated with the new conformation.
Returns the number of rotations performed.

# Examples
```julia-repl
julia> Mutators.Dihedral.apply!(state, mutator)
```
See also: [`Common.rotate_dihedral!`](@ref Common)
"""
function apply!(state::Common.State, mutator::MutatorConfig)
    
    count::Int64 = 0
    for dihedral in mutator.dihedrals
        if rand() < mutator.p_mut
            Common.rotate_dihedral!(state.xyz, dihedral, mutator.angle_sampler())
            count += 1
        end # end if
    end # end for
    return count
end # end function

end # end module