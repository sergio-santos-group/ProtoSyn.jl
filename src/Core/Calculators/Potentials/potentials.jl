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

    if update_forces & !(f === nothing)
        if selection !== nothing
            for (i, atom) in enumerate(selection(pose; gather = true))
                pose.state.f[:, atom.index] += f[:, i]
            end
        else
            for atom_index in 1:pose.state.size
                pose.state.f[:, atom_index] += f[:, atom_index]
            end
        end
    end

    pose.state.e[:Total] = e
    return e, f
end

# ------------------------------------------------------------------------------

# * Available potential functions

include("flat_bottom_potential.jl")
include("coulomb_potential.jl")