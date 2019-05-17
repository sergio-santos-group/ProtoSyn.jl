@doc raw"""
    evaluate!(st::Common.State, topology::Vector{DistanceFBR}, do_forces::Bool = false)

Evaluate an array of [Restraints.DistanceFBR](@ref Forcefield) using the current [`Common.State`](@ref),
calculate and update state.energy according to the equations defined in each stage of the flat-bottomed restraint.
If `do_forces` flag is set to `true`, calculate and update `state.forces`.
Return the component energy value (kJ mol⁻¹).

# Examples
```julia-repl
julia> Forcefield.Restraints.evaluate!(distances, state)
0.500
```
"""
function evaluate!(st::Common.State, topology::Vector{DistanceFBR}, do_forces::Bool = false)
    # All distances are in nm

    eDistanceFBR = .0
    v12          = [.0, .0, .0]
    f            = [.0, .0, .0]
    forces       = st.forces
    for pair in topology
        d12 = .0
        @inbounds for i=1:3
            delta  = st.xyz[pair.a2, i] - st.xyz[pair.a1, i]
            v12[i] = delta
            d12   += delta ^ 2
        end # end for
        d12 = sqrt(d12)
        if d12 < pair.r1
            dr1           = pair.r1 - pair.r2
            fconst        = pair.c * dr1
            e1            = fconst * dr1
            dr            = d12 - pair.r1
            eDistanceFBR += fconst * dr + e1 * 0.5
        elseif d12 < pair.r2
            dr            = d12 - pair.r2
            fconst        = pair.c * dr
            eDistanceFBR += fconst * dr * 0.5
        elseif d12 < pair.r3
            continue
        elseif d12 < pair.r4
            dr            = d12 - pair.r3
            fconst        = pair.c * dr
            eDistanceFBR += fconst * dr * 0.5
        else
            dr2           = pair.r4 - pair.r3
            fconst        = (pair.c * dr2)
            e2            = fconst * dr2
            dr            = d12 - pair.r4
            eDistanceFBR += fconst * dr + e2 * 0.5
        end # end if
        if do_forces
            @inbounds for i=1:3
                f = v12[i] * fconst / d12
                forces[pair.a1, i] += f
                forces[pair.a2, i] -= f
            end # end for
        end # end if
    end # end for

    Common.set_energy_component!(st.energy, :contacts, eDistanceFBR)
    return eDistanceFBR
end # end function


@doc raw"""
    evaluate!(st::Common.State, topology::Vector{DihedralFBR}, do_forces::Bool = false)

Evaluate an array of [Restraints.DihedralFBR](@ref Forcefield) using the current [`Common.State`](@ref),
calculate and update state.energy according to the equations defined in each stage of the flat-bottomed restraint.
If `do_forces` flag is set to `true`, calculate and update `state.forces`.
Return the component energy value (kJ mol⁻¹).

# Examples
```julia-repl
julia> Forcefield.Restraints.evaluate!(dihedrals, state)
0.500
```
"""
function evaluate!(st::Common.State, topology::Vector{DihedralFBR}, do_forces::Bool = false)
    # All distances are in nm and angles in rad

    eDihedralFBR = .0
    v12          = [.0, .0, .0]
    v32          = [.0, .0, .0]
    v34          = [.0, .0, .0]
    f1           = [.0, .0, .0]
    f3           = [.0, .0, .0]
    f4           = [.0, .0, .0]
    m            = [.0, .0, .0]
    n            = [.0, .0, .0]
    coords       = st.xyz
    forces       = st.forces

    for dihedral in topology

        d32Sq  = .0
        d3432  = .0
        d1232  = .0

        @inbounds for i=1:3
            a1 = coords[dihedral.a1, i]
            a2 = coords[dihedral.a2, i]
            a3 = coords[dihedral.a3, i]
            a4 = coords[dihedral.a4, i]

            delta12 = a2 - a1
            v12[i]  = delta12

            delta32 = a2 - a3
            v32[i]  = delta32
            d32Sq  += delta32*delta32

            delta34 = a4 - a3
            v34[i]  = delta34

            d3432  += delta34 * delta32
            d1232  += delta12 * delta32
        end # end for

        m       = cross(v12, v32)
        n       = cross(v32, v34)
        d32     = sqrt(d32Sq)
        phi     = - atan(d32 * dot(v12, n), dot(m, n))

        if phi <= dihedral.r1
            dr1           = dihedral.r1 - dihedral.r2
            cdr1          = dihedral.c * dr1
            e1            = cdr1 * dr1 * 0.5
            dr            = phi - dihedral.r1
            eDihedralFBR += cdr1 * dr + e1
            dVdphi_x_d32  = cdr1 * d32
        elseif phi <= dihedral.r2
            dr            = phi - dihedral.r2
            cdr2          = dihedral.c * dr
            eDihedralFBR += cdr2 * dr * 0.5 
            dVdphi_x_d32  = cdr2 * d32
        elseif phi <= dihedral.r3
            continue
        elseif phi <= dihedral.r4
            dr            = phi - dihedral.r3
            cdr4          = dihedral.c * dr
            eDihedralFBR += cdr4 * dr * 0.5
            dVdphi_x_d32  = cdr4 * d32
        else
            dr2           = dihedral.r4 - dihedral.r3
            cdr5          = dihedral.c * dr2
            e2            = cdr5 * dr2 * 0.5
            dr            = phi - dihedral.r4
            eDihedralFBR += cdr5 * dr + e2
            dVdphi_x_d32  = cdr5 * d32
        end # end if
        if do_forces
            f1mm = -dVdphi_x_d32 / dot(m, m)
            f4nn = dVdphi_x_d32 / dot(n, n)
            f3_1 = (d3432/d32Sq - 1.0)
            f3_2 = (d1232/d32Sq)
            @inbounds for i=1:3
                f1 = m[i] * f1mm
                f4 = n[i] * f4nn
                f3 = f4 * f3_1 - f1 * f3_2

                forces[dihedral.a1, i] += f1
                forces[dihedral.a2, i] += (-f1 - f3 - f4)
                forces[dihedral.a3, i] += f3
                forces[dihedral.a4, i] += f4
            end # end for
        end # end if
    end # end for

    Common.set_energy_component!(st.energy, :dihedralFBR, eDihedralFBR)
    return eDihedralFBR
end # end function