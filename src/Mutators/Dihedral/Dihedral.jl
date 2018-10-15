module Dihedral

using ..Common
using ..Aux

# @doc raw"""
#     ConfigParameters(p_mut::Float64 = 0.1)

# Define the runtime parameters for Dihedral movements.

# # Arguments
# - `p_mut::Float64`: Probability of mutation of this dihedral (Default: 0.1).

# # Examples
# ```julia-repl
# julia> Mutators.Diehdral.ConfigParameters(0.2)
# ConfigParameters(p_mut=0.2)

# julia> Mutators.Diehdral.ConfigParameters()
# ConfigParameters(p_mut=0.1)
# ```
# """
# mutable struct ConfigParameters

#     p_mut::Float64
#     step_size::Float64

#     ConfigParameters(; p_mut::Float64 = 0.1, step_size::Float64 = 0.1) = new(p_mut, step_size)
# end
# Base.show(io::IO, b::ConfigParameters) = print(io, "ConfigParameters(p_mut=$(b.p_mut), step_size=$(b.step_size))")

# -----------------------------------------------------------------------------------------------------------

mutable struct DihedralMutator
    dihedrals::Vector{MutableDihedral}
    pmut::Float64
    angle_sampler::Function
    stepsize::Float64
end

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
    end
end

end