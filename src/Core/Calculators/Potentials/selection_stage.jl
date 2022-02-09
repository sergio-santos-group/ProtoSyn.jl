# * No selection

# Main function
function resolve_selection(A::Type{<: ProtoSyn.AbstractAccelerationType},
    pose::Pose,
    potential::Function,
    update_forces::Bool,
    verlet_list::Union{VerletList, Nothing},
    selection::Nothing,
    mask::MaskMap)

    coords = pose.state.x.coords[:]
    return resolve_calculation(A, pose, potential, update_forces, verlet_list, coords, mask, collect(1:size(pose.state.x.coords)[2]))
end

# ! ---

# * With Selection

# Main function
function resolve_selection(A::Type{<: ProtoSyn.AbstractAccelerationType},
    pose::Pose,
    potential::Function,
    update_forces::Bool,
    verlet_list::Union{VerletList, Nothing},
    selection::AbstractSelection,
    mask::Opt{Union{ProtoSyn.Mask{<: ProtoSyn.AbstractContainer}, Matrix{<: AbstractFloat}}})

    sele = ProtoSyn.promote(selection, Atom)(pose).content
    N = count(sele)
    if N == 0
        @warn "The provided selection resulted in 0 atoms selected to apply potential."
        return 0.0, nothing
    end
    coords = pose.state.x.coords[:, sele][:]
    indexes = findall(==(true), sele)

    return resolve_calculation(A, pose, potential, update_forces, verlet_list, coords, mask, indexes)
end

function resolve_selection(A::Type{<: ProtoSyn.AbstractAccelerationType},
    pose::Pose,
    potential::Function,
    update_forces::Bool,
    verlet_list::Union{VerletList, Nothing},
    selection::AbstractSelection,
    mask::Function)

    sele = ProtoSyn.promote(selection, Atom)(pose).content
    N = count(sele)
    if N == 0
        @warn "The provided selection resulted in 0 atoms selected to apply potential."
        return 0.0, nothing
    end
    coords = pose.state.x.coords[:, sele][:]
    indexes = findall(==(true), sele)

    return resolve_calculation(A, pose, potential, update_forces, verlet_list, coords, mask(pose, selection), indexes)
end