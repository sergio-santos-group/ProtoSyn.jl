using Statistics
using LinearAlgebra

export align!

"""
    align!(mobile::Pose, target::Pose)
    align!(mobile::Pose, target::Pose, selection::ProtoSyn.AbstractSelection)
    align!(mobile::Pose, target::Pose, mobile_selection::ProtoSyn.AbstractSelection, target_selection::ProtoSyn.AbstractSelection)
    
Application of the Kabsch algorithm. Applies a rotation + translation movement
on the `mobile` [`Pose`](@ref) instance in order to align to the `target`
[`Pose`](@ref) instance. If a `selection` is provided, only the subset of
selected [`Atom`](@ref) instances (on both [`Pose`](@ref) structures) will be
considered to calculate the necessary rotation + translation movement
(minimizing the RMSD). If two `AbstractSelection` instances are provided
(`mobile_selection` and `target_selection`), each is applied to the respective
[`Pose`](@ref) instances (`mobile` and `target`, respectively) in order to
calculate the necessary rotation + translation movement. Sets `mobile.state.c2i`
to `true` and returns the altered `mobile` [`Pose`](@ref) instance.

!!! ukw "Note:"
    This function can also align [`Fragment`](@ref) instances.

# See also
[`rmsd`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.align!(pose, pose_mod)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 343
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)

julia> ProtoSyn.align!(pose, pose_mod, an"CA")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 343
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)

julia> ProtoSyn.align!(pose, pose_mod, an"CA", an"CB")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 343
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function align!(mobile::Pose, target::Pose, selection::AbstractSelection)
    
    sele            = promote(selection, Atom)
    mobile_coords   = copy(mobile.state.x.coords[:, sele(mobile).content])
    target_coords   = copy(target.state.x.coords[:, sele(target).content])

    c_mobile                 = mean(mobile_coords, dims = 2)
    c_target                 = mean(target_coords, dims = 2)
    _mobile_coords           = mobile_coords .- c_mobile
    _target_coords           = target_coords .- c_target
       
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


"""
    center_of_mass(pose::Pose)
    
Return the center of mass X, Y and Z cartesian coordinates of the given
[`Pose`](@ref) `pose` (based on the current cartesian coordinates - make sure
the [`Pose`](@ref) `pose` is synched, using the [`sync!`](@ref) method).

    center_of_mass(pose::Pose, selection::AbstractSelection)

Return the center of mass X, Y and Z cartesian coordinates of the given
[`Pose`](@ref) `pose`, taking into consideration only the subset of selected
[`Atom`](@ref) instances in the `AbstractSelection` `selection` (based on the
current cartesian coordinates - make sure the [`Pose`](@ref) `pose` is synched,
using the [`sync!`](@ref) method).

    center_of_mass(pose::Pose, idxs::Vector{Int})

Return the center of mass X, Y and Z cartesian coordinates of the given
[`Pose`](@ref) `pose`, taking into consideration only the subset of 
[`Atom`](@ref) instances in the vector `idxs` (by [`Atom`](@ref) index, based on
the current cartesian coordinates - make sure the [`Pose`](@ref) `pose` is
synched, using the [`sync!`](@ref) method).

# Examples
```
julia> ProtoSyn.center_of_mass(pose)
3×1 Matrix{Float64}:
 39.85855147920603
 14.995282315671613
 -0.016516024315774907

julia> ProtoSyn.center_of_mass(pose, an"CA")
3×1 Matrix{Float64}:
 37.56949530961898
 14.249844760318357
 -5.4078476622375185e-16
```
"""
function center_of_mass(pose::Pose)
    return mean(pose.state.x.coords, dims = 2)
end

function center_of_mass(pose::Pose, selection::AbstractSelection)
    idxs = findall(promote(selection, Atom)(pose).content)
    return center_of_mass(pose, idxs)
end

center_of_mass(pose::Pose, selection::Nothing) = center_of_mass(pose)

function center_of_mass(pose::Pose, idxs::Vector{Int})
    return mean(pose.state.x[:, idxs], dims = 2)
end