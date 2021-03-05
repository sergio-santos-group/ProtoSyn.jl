using Statistics
using LinearAlgebra

export align!

"""
    align!(mobile::Pose, target::Pose)
    align!(mobile::Pose, target::Pose, selection::ProtoSyn.AbstractSelection)
    align!(mobile::Pose, target::Pose, mobile_selection::ProtoSyn.AbstractSelection, target_selection::ProtoSyn.AbstractSelection)
    
Application of the Kabsch algorithm. Applies a rotation + translation movement
on the mobile Pose in order to align to the target Pose. If a selection is
provided, only the subset of selected Atoms (on both Poses) will be considered
to calculate the necessary rotation + translation movement (minimizing the
RMSD). If two selection are provided, each is applied to the respective Pose.
Sets `mobile.state.c2i` to `true`.
!!! note
    This function can align Fragments.

# Examples
```jldoctest
julia> ProtoSyn.align!(pose2, pose1)

julia> ProtoSyn.align!(pose2, pose1, an"CA")

julia> ProtoSyn.align!(pose2, pose1, an"CA", an"CB")
```
"""
function align!(mobile::Pose, target::Pose, selection::AbstractSelection)
    
    sele            = promote(selection, Atom)
    mobile_coords   = copy(mobile.state.x.coords[:, sele(mobile).content])
    target_coords   = copy(target.state.x.coords[:, sele(target).content])

    c_mobile                 = mean(mobile_coords, dims = 2)
    c_target                 = mean(target_coords, dims = 2)
    _mobile_coords         .-= c_mobile
    _target_coords         .-= c_target
       
    cm                       = _mobile_coords * _target_coords'
    u, d, vt                 = svd(cm)
    rot_matrix               = vt * u'
       
    centered_mobile          = mobile.state.x.coords .- c_mobile
    mobile.state.x[:, 1:end] = (centered_mobile' * rot_matrix')' .+ c_target

    return mobile
end


function align!(mobile::Pose, target::Pose, mobile_selection::AbstractSelection, target_selection::AbstractSelection)
    
    # 1) Apply selection
    _mobile_selection        = promote(mobile_selection, Atom)
    mobile_mask              = _mobile_selection(mobile).content
    mobile_coords            = copy(mobile.state.x.coords[:, mobile_mask])
       
    _target_selection        = promote(target_selection, Atom)
    target_mask              = _target_selection(target).content
    target_coords            = copy(target.state.x.coords[:, target_mask])
       
    c_mobile                 = mean(mobile_coords, dims = 2)
    c_target                 = mean(target_coords, dims = 2)
    mobile_coords          .-= c_mobile
    target_coords          .-= c_target
       
    cm                       = mobile_coords * target_coords'
    u, d, vt                 = svd(cm)
    rot_matrix               = vt * u'
       
    centered_mobile          = mobile.state.x.coords .- c_mobile
    mobile.state.x[:, 1:end] = (centered_mobile' * rot_matrix')' .+ c_target

    return mobile
end


align!(mobile::Pose, target::Pose) = begin
    return align!(mobile, target, ProtoSyn.TrueSelection{Atom}())
end