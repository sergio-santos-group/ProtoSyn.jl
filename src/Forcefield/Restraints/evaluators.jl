#TODO: Document function
# function calc_eSol!(st::Common.State, residues::Vector{Common.Residue})
    
#     # Create solvation pairs
#     d = Dict("Q" => -3.5, "W" => -0.9, "T" => -0.7, "C" =>  2.5, "P" => -1.6, "V" =>  4.2, "L" =>  3.8, "M" =>  1.9, "N" => -3.5, "H" => -3.2,
#              "A" =>  1.8, "D" => -3.5, "G" => -0.4, "E" => -3.5, "Y" => -1.3, "I" =>  4.5, "S" => -0.8, "K" => -3.9, "R" => -4.5, "F" =>  2.8)
#     solv_pairs = map(x -> SolvPair(x.cα, d[Aux.conv321(string(x.name))]), residues)
    
#     # Calculate solvation energy
#     n_atoms = length(solv_pairs)
#     l = zeros(Float64, 3)
#     e_sol::Float64 = 0.0
#     sum_f::Float64 = 0.0
#     for i in 1:(n_atoms - 1)
#         Ai = @view st.xyz[solv_pairs[i].i,:]
#         sum_f = 0.0
#         for j in (i+1):n_atoms
#             Aj = @view st.xyz[solv_pairs[j].i,:]
#             @. l[:] = Aj - Ai
#             dIJ = norm(l)
#             f = 1.0 / (1.0 + exp(-(12.0-dIJ) / 0.4))
#             sum_f += f
#         end
#         if ((sum_f < 21.0) && (solv_pairs[i].coef > 0.0)) || ((sum_f > 21.0) && (solv_pairs[i].coef < 0.0))
#             e_sol += solv_pairs[i].coef * (21.0 - sum_f)
#         end
#     end

#     st.energy.comp["eSol"] = e_sol
#     st.energy.eTotal = e_sol
#     return e_sol
# end


@doc raw"""
    evaluate!(topology::Vector{DistanceFBR}, state::Common.State[, do_forces::Bool = false])::Float64

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
function evaluate!(topology::Vector{DistanceFBR}, st::Common.State; do_forces::Bool = false)
    # All distances are in nm

    eDistanceFBR::Float64 = 0.0
    v12::Vector{Float64} = zeros(Float64, 3)
    f::Vector{Float64} = zeros(Float64, 3)
    for pair in topology
        @views @. v12 = st.xyz[pair.a2, :] - st.xyz[pair.a1, :]
        d12 = norm(v12)
        if d12 < pair.r1
            dr1::Float64 = pair.r1 - pair.r2
            e1::Float64 = pair.c * dr1 * dr1 * 0.5
            dr = d12 - pair.r1
            eDistanceFBR += (pair.c * dr1) * dr + e1
            fconst = (pair.c * dr1)
        elseif d12 < pair.r2
            dr = d12 - pair.r2
            eDistanceFBR += pair.c * dr * dr * 0.5
            fconst = (pair.c * dr)
        elseif d12 < pair.r3
            continue
        elseif d12 < pair.r4
            dr = d12 - pair.r3
            eDistanceFBR += pair.c * dr * dr * 0.5
            fconst = (pair.c * dr)
        else
            dr2::Float64 = pair.r4 - pair.r3
            e2::Float64 = pair.c * dr2 * dr2 * 0.5
            dr = d12 - pair.r4
            eDistanceFBR += (pair.c * dr2) * dr + e2
            fconst = (pair.c * dr2)
        end
        if do_forces
            @. f = v12 * fconst / d12
            @. st.forces[pair.a1, :] += f
            @. st.forces[pair.a2, :] -= f
        end
    end

    st.energy.comp["eDistanceFBR"] = eDistanceFBR
    st.energy.eTotal = eDistanceFBR
    return eDistanceFBR
end


@doc raw"""
    evaluate!(topology::Vector{DihedralFBR}, state::Common.State[, do_forces::Bool = false])::Float64

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
function evaluate!(topology::Vector{DihedralFBR}, st::Common.State; do_forces::Bool = false)
    # All distances are in nm and angles in rad

    eDihedralFBR::Float64 = 0.0
    v12 = zeros(Float64, 3)
    v32 = zeros(Float64, 3)
    v34 = zeros(Float64, 3)
    f1 = zeros(Float64, 3)
    f3 = zeros(Float64, 3)
    f4 = zeros(Float64, 3)
    m = zeros(Float64, 3)
    n = zeros(Float64, 3)

    for dihedral in topology
        # println(dihedral)
        @views @. v12 = st.xyz[dihedral.a2, :] - st.xyz[dihedral.a1, :]
        @views @. v32 = st.xyz[dihedral.a2, :] - st.xyz[dihedral.a3, :]
        @views @. v34 = st.xyz[dihedral.a4, :] - st.xyz[dihedral.a3, :]
        m = cross(v12, v32)
        n = cross(v32, v34)
        d32Sq = dot(v32, v32)
        d32 = sqrt(d32Sq)
        phi = - atan(d32 * dot(v12, n), dot(m, n))
        # println("Current: $(rad2deg(phi))")

        if phi <= dihedral.r1
            dr1::Float64 = dihedral.r1 - dihedral.r2
            e1::Float64 = dihedral.c * dr1 * dr1 * 0.5
            dr = phi - dihedral.r1
            eDihedralFBR += (dihedral.c * dr1) * dr + e1
            dVdphi_x_d32 = (dihedral.c * dr1) * d32
        elseif phi <= dihedral.r2
            dr = phi - dihedral.r2
            eDihedralFBR += dihedral.c * dr * dr * 0.5 
            dVdphi_x_d32 = dihedral.c * dr * d32
        elseif phi <= dihedral.r3
            continue
        elseif phi <= dihedral.r4
            dr = phi - dihedral.r3
            eDihedralFBR += dihedral.c * dr * dr * 0.5
            dVdphi_x_d32 = dihedral.c * dr * d32
        else
            dr2::Float64 = dihedral.r4 - dihedral.r3
            e2::Float64 = dihedral.c * dr2 * dr2 * 0.5
            dr = phi - dihedral.r4
            eDihedralFBR += (dihedral.c * dr2) * dr + e2
            dVdphi_x_d32 = (dihedral.c * dr2) * d32
        end
        if do_forces
            f1 .= m .* (-dVdphi_x_d32 / dot(m, m))
            f4 .= n .* ( dVdphi_x_d32 / dot(n, n))
            f3 .= f4 .* (dot(v34, v32)/d32Sq - 1.0) .- f1 .* (dot(v12, v32)/d32Sq)
            # println("Forces:\n a1: $(f1)\n a1: $(-f1-f3-f4)\n a1: $(f3)\n a1: $(f4)")
            @. st.forces[dihedral.a1, :] += f1
            @. st.forces[dihedral.a2, :] += (-f1 - f3 - f4)
            @. st.forces[dihedral.a3, :] += f3
            @. st.forces[dihedral.a4, :] += f4
        end
    end

    st.energy.comp["eDihedralFBR"] = eDihedralFBR
    st.energy.eTotal = eDihedralFBR
    return eDihedralFBR
end