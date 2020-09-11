using LinearAlgebra

export distance

"""
    distance(at1::AtomState, at2::AtomState)
    
Calculates the distance between the two `AtomState` instances, based on the
cartesian coordinates. Note: Make sure the corresponding pose has been synced.
Returns result in ???.

# Examples
```jldoctest
julia> d = distance(pose.state[1], pose.state[2])
1.004
```
"""
function distance(at1::AtomState, at2::AtomState)
    return norm(at1.t - at2.t)
end

export angle
function angle(at1::AtomState, at2::AtomState, at3::AtomState)

    v21 = at1.t - at2.t
    v23 = at3.t - at2.t
    return acos(dot(v21, v23) / (norm(v21) * norm(v23)))
end


export dihedral
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