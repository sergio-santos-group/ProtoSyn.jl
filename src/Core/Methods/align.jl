using Statistics
using LinearAlgebra

export align!

"""
    align!(mobile::Pose, target::Pose)
    align!(mobile::Pose, target::Pose, selection::ProtoSyn.AbstractSelection)
    
Application of the Kabsch algorithm. Applies a rotation + translation movement
on the mobile Pose in order to align to the target Pose. If a selection is
provided, only the subset of selected Atoms will be considered to calculate the
necessary rotation + translation movement (minimizing the RMSD). Sets
`mobile.state.c2i` to `true`.

# Examples
```jldoctest
julia> ProtoSyn.align!(pose2, pose1)

julia> ProtoSyn.align!(pose2, pose1, an"CA")
```
"""
function align!(mobile::Pose, target::Pose, selection::AbstractSelection)
    
    # 1) Apply selection
    sele            = promote(selection, Atom)
    mobile_coords   = copy(mobile.state.x.coords[:, sele(mobile).content])
    target_coords   = copy(target.state.x.coords[:, sele(target).content])

    # 2) Get centroid and rotation matrix
    # 2.1) Center each centroid on [0, 0, 0]
    c_mobile        = mean(mobile_coords, dims = 2)
    c_target        = mean(target_coords, dims = 2)
    mobile_coords .-= c_mobile
    target_coords .-= c_target

    # 2.2) Obtain correlation matrix
    cm              = mobile_coords * target_coords'
    u, d, vt        = svd(cm)
    rot_matrix      = vt * u'

    # 3) Apply transformation (rotation+ translation)
    centered_mobile = mobile.state.x.coords .- c_mobile
    mobile.state.x[:, 1:end] = (centered_mobile' * rot_matrix')' .+ c_target

    return mobile
end

align!(mobile::Pose, target::Pose) = begin
    return align!(mobile, target, ProtoSyn.TrueSelection{Atom}())
end