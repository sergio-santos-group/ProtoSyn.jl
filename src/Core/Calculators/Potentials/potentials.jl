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

apply_potential!(A::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, potential::Function, update_forces::Bool, verlet_list::Union{VerletList, Nothing}, selection::Union{Nothing, ProtoSyn.AbstractSelection}, mask::MaskMap) = begin
    e, f = apply_potential(A, pose, potential, update_forces, verlet_list, selection, mask)

    return e, f
end

# ------------------------------------------------------------------------------

# * Available potential functions

include("flat_bottom_potential.jl")
include("coulomb_potential.jl")