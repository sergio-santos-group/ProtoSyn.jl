module Dihedral

using ..Common
using ..Aux
using ..Mutators

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
@Base.kwdef mutable struct DihedralMutator <: Mutators.AbstractMutator
    dihedrals::Vector{Common.Dihedral}
    angle_sampler::Function
    p_mut::Float64 = 0.0
    step_size::Float64 = 0.0
    # ac_ratio::Float64 = 0.0
end
# DihedralMutator(dihedrals::Vector{Common.Dihedral}, angle_sampler::Function; p_mut = 0.0, step_size = 0.0) = DihedralMutator(dihedrals, angle_sampler, p_mut, step_size)
Base.show(io::IO, b::DihedralMutator) = print(io, "DihedralMutator(dihedrals=$(length(b.dihedrals)), angle_sampler=$(string(b.angle_sampler)), p_mut=$(b.p_mut), step_size=$(b.step_size))")


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
function apply!(state::Common.State, mutator::DihedralMutator)
    
    count::Int64 = 0
    for dihedral in mutator.dihedrals
        if rand() < mutator.p_mut
            Common.rotate_dihedral!(state.xyz, dihedral, mutator.angle_sampler())
            count += 1
        end
    end
    # Adjust step_size: 0.2 is the desired acceptance ratio, 0.001 is the minimum step_size <- Are these modifiable ?
    # δ = mutator.ac_ratio > 0.2 ? 1.05 : 0.5
    # mutator.step_size = max(0.001, min(mutator.step_size * δ, π))
    return count
end

end

# adjust_step_size_rfnm = @Common.callback 1 function _adjust_step_size_refnm(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.Driver, ac::Float64, args...)
#     d::Float64 = 1.0
#     if ac > acceptance_ratio + ar_buffer_zone
#         d = 1.05
#     elseif ac < acceptance_ratio - ar_buffer_zone
#         d = 0.95
#     end
#     # soft_dh_mutator.step_size = max(min_step_s, min(soft_dh_mutator.step_size * d, π))
#     # soft_cs_mutator.step_size = max(min_step_s, min(soft_cs_mutator.step_size * d, π))
#     soft_br_mutator.step_size = max(min_step_s, min(soft_br_mutator.step_size * d, π))
# end