module Dihedral

using ..Common
using ..Aux

#TODO: Document structure
mutable struct DihedralMutator
    dihedrals::Vector{Common.Dihedral}
    pmut::Float64
    angle_sampler::Function
    stepsize::Float64
end

#TODO: Document function
@doc raw"""
    run!(state::Common.State, dihedrals::Array{NewDihedral, 1}, params::Config.Parameters, angle_sampler::Function[, ostream::IO = stdout])

Iterate over a list of [`NewDihedral`](@ref) (`dihedrals`) and perform dihedral movements on the current
[`Common.State`](@ref). The probability of each dihedral undergo movements is defined in the
[`ConfigParameters`](@ref) (`params.pmut`). The new angle is obtained from `angle_sampler`, who should
return a `Float64` in radians.
After movement, the [`Common.State`](@ref) is updated with the new conformation.
Any logging is written to `ostream` (Default: `stdout`).

# Examples
```julia-repl
julia> Mutators.Diehdral.run!(state, dihedrals, params, () -> randn())
```
See also: [`rotate_dihedral!`](@ref)
"""
@inline function run!(state::Common.State, mutator::DihedralMutator)
    for dihedral in mutator.dihedrals
        if rand() < mutator.pmut
            Common.rotate_dihedral!(state.xyz, dihedral, mutator.angle_sampler())
        end
        xyz[idxs, :] = (rmat * (xyz[idxs, :] .- pivot)')' .+ pivot
    end
end





end