using LinearAlgebra

export distance
"""
    distance(at1::AtomState, at2::AtomState)
    
Calculates the distance between the two [`AtomState`](@ref) instances (`at1` and
`at2`), based on the cartesian coordinates. Note: Make sure the corresponding
[`Pose`](@ref) instance has been synched (using the [`sync!`](@ref) method).
Returns result in Angstrom (Å).

# See also
[`angle`](@ref) [`dihedral`](@ref)

# Examples
```jldoctest
julia> d = ProtoSyn.distance(pose.state[1], pose.state[2])
1.0093879999999997
```
"""
function distance(at1::AtomState, at2::AtomState)
    return norm(at1.t - at2.t)
end

export angle
"""
    angle(at1::AtomState, at2::AtomState, at3::AtomState)
    
Calculates the angle between the three [`AtomState`](@ref) instances (`at1`,
`at2` and `at3`), based on the cartesian coordinates. Note: Make sure the
corresponding [`Pose`](@ref) instance has been synched (using the
[`sync!`](@ref) method). Returns result in radians.

# See also
[`distance`](@ref) [`dihedral`](@ref)

# Examples
```jldoctest
julia> a = ProtoSyn.angle(pose.state[1], pose.state[2], pose.state[3])
0.6444967663659441
```
"""
function angle(at1::AtomState, at2::AtomState, at3::AtomState)

    v21 = at1.t - at2.t
    v23 = at3.t - at2.t
    a = dot(v21, v23) / (norm(v21) * norm(v23))
    return acos(ceil(a, digits = 15)) # ceil prevents erros when ∠ ≈ 180°
end


export dihedral
"""
    dihedral(at1::AtomState, at2::AtomState, at3::AtomState, at4::AtomState)
    
Calculates the dihedral angle between the four [`AtomState`](@ref) instances
(`at1`, `at2`, `at3` and `at4`), based on the cartesian coordinates. Note: Make
sure the corresponding [`Pose`](@ref) instance has been synched (using the
[`sync!`](@ref) method). Returns result in radians.

# Examples
```jldoctest
julia> a = ProtoSyn.dihedral(pose.state[1], pose.state[2], pose.state[3], pose.state[4])
-1.2318251145557122
```
"""
function dihedral(at1::AtomState, at2::AtomState, at3::AtomState, at4::AtomState)

    b1 = at2.t - at1.t
    b2 = at3.t - at2.t
    b3 = at4.t - at3.t
    n1 = cross(b1, b2)
    n2 = cross(b2, b3)
    x = dot(cross(n1, n2), b2) / sqrt(dot(b2, b2))
    y = dot(n1, n2)
    return atan(x, y) # in rad
end


export rmsd
"""
    rmsd(pose1::Pose, pose2::Pose)
    rmsd(pose1::Pose, pose2::Pose, selection::AbstractSelection)
    
Calculates the RMSD value between 2 [`Pose`](@ref) instances, based on the
cartesian coordinates. Note: Make sure the poses have been synched beforehand
(using the [`sync!`](@ref) method). If an `AbstractSelection` `selection` is
provided, calculate the RMSD values of only the selected subset of
[`Atom`](@ref) instances. Returns RMSD result in Angstrom (Å). 

# See also
[`align!`](@ref) [`getdihedral`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.rmsd(pose, pose_mod)
34.443215682826676

julia> ProtoSyn.rmsd(pose, pose_mod, an"CA")
32.5063913965703
```
"""
function rmsd(pose1::Pose, pose2::Pose, selection::AbstractSelection)
    sele = promote(selection, Atom)
    p1 = pose1.state.x.coords[:, sele(pose1).content]
    p2 = pose2.state.x.coords[:, sele(pose2).content]
    d = p1 .- p2
    return sqrt(sum(d.*d)/size(d)[2])
end


function rmsd(pose1::Pose, pose2::Pose)
    d = pose1.state.x.coords .- pose2.state.x.coords
    return sqrt(sum(d.*d)/size(d)[2])
end