module Dihedral

using ..Common
using ..Aux

@doc raw"""
    DihedralMutator(dihedrals::Vector{Common.Dihedral}, p_mut::Float64, angle_sampler::Function, step_size::Float64)

Holds all necessary parameters for the correct simulation of dihedral movements.

# Arguments
- `dihedrals::Vector{Common.Dihedral}`: List of dihedrals avaliable to be rotated in dihedral movements.
- `p_mut::Float64`: Probability of rotation of each dihedral.
- `angle_sampler::Function`: Function responsible for defining the new angle for the dihedral. Should return a `Float64`.
- `step_size::Float64`: Scalar that defines the amount of change resulting for a dihedral movement.

# Examples
```julia-repl
julia> Mutators.Diehdral.DihedralMutator(dihedrals, 0.05, randn, 0.25)
DihedralMutator(dihedrals=68, p_pmut=0.05, angle_sampler="randn", step_size=0.25)
```
See also: [`run!`](@ref)
"""
mutable struct DihedralMutator
    dihedrals::Vector{Common.Dihedral}
    p_mut::Float64
    angle_sampler::Function
    step_size::Float64
end
Base.show(io::IO, b::DihedralMutator) = print(io, "DihedralMutator(dihedrals=$(length(b.dihedrals)), p_mut=$(b.p_mut), angle_sampler=$(string(b.angle_sampler)), step_size=$(b.step_size))")

#TODO: Document function
@doc raw"""
    run!(state::Common.State, dihedrals::Array{NewDihedral, 1}, params::Config.Parameters, angle_sampler::Function[, ostream::IO = stdout])

Iterate over a list of [`Common.Dihedral`](@ref) (`dihedrals`) and perform dihedral movements on the current
[`Common.State`](@ref). The probability of each dihedral undergo movements is defined in the
[`DihedralMutator`](@ref).`p_mut`. The new angle is obtained from [`DihedralMutator`](@ref).`angle_sampler`, who should
return a `Float64` in radians.
After movement, the [`Common.State`](@ref) is updated with the new conformation.

# Examples
```julia-repl
julia> Mutators.Diehdral.run!(state, mutator)
```
See also: [`Common.rotate_dihedral!`](@ref Common)
"""
@inline function run!(state::Common.State, mutator::DihedralMutator)
    for dihedral in mutator.dihedrals
        if rand() < mutator.p_mut
            Common.rotate_dihedral!(state.xyz, dihedral, mutator.angle_sampler())
        end
        # xyz[idxs, :] = (rmat * (xyz[idxs, :] .- pivot)')' .+ pivot
    end
end

end