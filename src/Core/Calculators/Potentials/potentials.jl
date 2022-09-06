using LinearAlgebra

# * Definition of useful alias

MaskMap   = Opt{Union{ProtoSyn.Mask{<: ProtoSyn.AbstractContainer}, Matrix{<: AbstractFloat}, Function}}

include("mask_stage.jl")
include("calculation_stage_sisd.jl")
include("calculation_stage_simd.jl")
include("calculation_stage_cuda.jl")
include("selection_stage.jl")

# ------------------------------------------------------------------------------

# * Variant dispatch

apply_potential(A::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, potential::Function, update_forces::Bool, verlet_list::Union{VerletList, Nothing}, selection::Union{Nothing, ProtoSyn.AbstractSelection}, mask::MaskMap) = resolve_selection(A, pose, potential, update_forces, verlet_list, selection, mask)
apply_potential(A::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, potential::Function, update_forces::Bool, verlet_list::Union{VerletList, Nothing}, mask::MaskMap) = resolve_selection(A, pose, potential, update_forces, verlet_list, nothing, mask)
apply_potential(A::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, potential::Function, update_forces::Bool, verlet_list::Union{VerletList, Nothing}, selection::Union{Nothing, ProtoSyn.AbstractSelection}) = resolve_selection(A, pose, potential, update_forces, verlet_list, selection, nothing)
apply_potential(A::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, potential::Function, update_forces::Bool, verlet_list::Union{VerletList, Nothing}) = resolve_selection(A, pose, potential, update_forces, verlet_list, nothing, nothing)
apply_potential(pose::Pose, potential::Function, update_forces::Bool, verlet_list::Union{VerletList, Nothing}, selection::Union{Nothing, ProtoSyn.AbstractSelection}, mask::MaskMap) = resolve_selection(ProtoSyn.acceleration.active, pose, potential, update_forces, verlet_list, selection, mask)
apply_potential(pose::Pose, potential::Function, update_forces::Bool, verlet_list::Union{VerletList, Nothing}, mask::MaskMap) = resolve_selection(ProtoSyn.acceleration.active, pose, potential, update_forces, verlet_list, nothing, mask)
apply_potential(pose::Pose, potential::Function, update_forces::Bool, verlet_list::Union{VerletList, Nothing}, selection::Union{Nothing, ProtoSyn.AbstractSelection}) = resolve_selection(ProtoSyn.acceleration.active, pose, potential, update_forces, verlet_list, selection, nothing)
apply_potential(pose::Pose, potential::Function, update_forces::Bool, verlet_list::Union{VerletList, Nothing}) = resolve_selection(ProtoSyn.acceleration.active, pose, potential, update_forces, verlet_list, nothing, nothing)


"""
    apply_potential!([A::Type{<: ProtoSyn.AbstractAccelerationType}], pose::Pose, potential::Function, update_forces::Bool, [verlet_list::Union{VerletList, Nothing}], [selection::Union{Nothing, ProtoSyn.AbstractSelection}], [mask::MaskMap])

Apply the given `potential` to the provided [`Pose`](@ref) `pose`. If the
`update_forces` flag is set to true, also calculate the forces acting on the
system. The call to this function can be further customized by providing an
optional [`VerletList`](@ref), an `AbstractSelection` `selection` or a `MaskMap`
`map` (this can be a [`Mask`](@ref), a `Matrix` or a `Function`).

# Examples
```
julia> ProtoSyn.Calculators.apply_potential!(ProtoSyn.SISD_0, pose, (d; kwargs...) -> 0.0, false, nothing, nothing, nothing)
(0.0, [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0])
```
"""
apply_potential!(A::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, potential::Function, update_forces::Bool, verlet_list::Union{VerletList, Nothing}, selection::Union{Nothing, ProtoSyn.AbstractSelection}, mask::MaskMap) = begin
    e, f = apply_potential(A, pose, potential, update_forces, verlet_list, selection, mask)

    return e, f
end

# ------------------------------------------------------------------------------

# * Available potential functions

include("flat_bottom_potential.jl")
include("coulomb_potential.jl")
include("bump_potential.jl")
include("harmonic_potential.jl")
include("lennard_jones_potential.jl")