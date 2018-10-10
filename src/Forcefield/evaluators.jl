@doc raw"""
    evaluate!(bonds::Array{Forcefield.HarmonicBond}, state::Common.State[, do_forces::Bool = false])::Float64
"""
function evaluate!(bonds::Array{Forcefield.HarmonicBond}, state::Common.State;
    do_forces::Bool = false)

    energy = 0.0
    v12 = zeros(Float64, 3)

    for bond in bonds
        @views @. v12 = state.xyz[bond.a2, :] - state.xyz[bond.a1, :]
        d12 = norm(v12)
        dr = d12 - bond.b0
        energy += bond.k * dr * dr
        if do_forces
            @. v12 *= (bond.k * dr / d12)
            @. state.forces[bond.a1, :] += v12
            @. state.forces[bond.a2, :] -= v12
        end
    end
    
    state.energy.eBond = 0.5 * energy
    
    return 0.5 * energy
end

@doc raw"""
    evaluate!(angles::Array{Forcefield.HarmonicAngle}, state::Common.State, do_forces::Bool = false)::Float64
"""
function evaluate!(angles::Array{Forcefield.HarmonicAngle}, state::Common.State;
    do_forces::Bool = false)

    v12 = zeros(Float64, 3)
    v32 = zeros(Float64, 3)
    f1 = zeros(Float64, 3)
    f3 = zeros(Float64, 3)
    energy = 0.0

    for angle in angles
        @views @. v12 = state.xyz[angle.a2, :] - state.xyz[angle.a1, :]
        @views @. v32 = state.xyz[angle.a2, :] - state.xyz[angle.a3, :]
        
        d12Sq = dot(v12, v12)
        d32Sq = dot(v32, v32)
        d12xd32 = sqrt(d12Sq * d32Sq)
        ctheta = dot(v12, v32) / d12xd32
        dtheta = acos(ctheta) - angle.θ

        energy += angle.k * dtheta * dtheta

        if do_forces
            fc = angle.k * dtheta / sqrt(1.0 - ctheta * ctheta)
            @. f1 = fc * (ctheta * v12/d12Sq - v32/d12xd32)
            @. f3 = fc * (ctheta * v32/d32Sq - v12/d12xd32)
            @. state.forces[angle.a1, :] += f1
            @. state.forces[angle.a3, :] += f3
            @. state.forces[angle.a2, :] -= (f1 + f3)
        end
    end

    state.energy.eAngle = 0.5 * energy
    
    return 0.5 * energy
end


# export evaldihedralRB
# function evaldihedralRB(dihedralsRB::Array{FFComponents.DihedralRB}, xyz::Array{Float64, 2}, forces::Array{Float64, 2};
#     do_forces = false)

#     energy = 0.0
#     for dihedral in dihedralsRB
#         v12 = xyz[dihedral.a2, :] - xyz[dihedral.a1, :]
#         v32 = xyz[dihedral.a2, :] - xyz[dihedral.a3, :]
#         v34 = xyz[dihedral.a4, :] - xyz[dihedral.a3, :]
#         m = cross(v12, v32)
#         n = cross(v32, v34)
#         d32Sq = norm(v32) ^ 2
#         d32 = norm(v32)
#         phi = atan2(d32 * dot(v12, n), dot(m, n))
#         cpsi = cos(phi - pi)
        
#         energy += dihedral.c0 + cpsi * (dihedral.c1 + cpsi * (dihedral.c2 + cpsi * (dihedral.c3 + cpsi * (dihedral.c4 + cpsi * dihedral.c5))))
        
#         if do_forces
#             cphi = -cpsi
#             dVdphi_x_d32 = sin(phi) * d32 * (dihedral.c1 + cphi * (-2.0 * dihedral.c2 + cphi * (3.0 * dihedral.c3 + cphi * (-4.0 * dihedral.c4 + cphi * 5.0 * dihedral.c5))))
#             f1 = m .* -dVdphi_x_d32/dot(m, m)
#             f4 = n .*  dVdphi_x_d32/dot(n, n)
#             f3 = -f4
#             f3 -= f1 .* dot(v12, v32)/d32Sq
#             f3 += f4 .* dot(v34, v32)/d32Sq
#             forces[dihedral.a1, :] += f1
#             forces[dihedral.a2, :] += -f1 - f3 - f4
#             forces[dihedral.a3, :] += f3
#             forces[dihedral.a4, :] += f4
#         end
#     end
#     return energy
# end

@doc raw"""
    evaluate!(dihedralsCos::Array{Forcefield.DihedralCos}, state::Common.State, do_forces::Bool = false)::Float64
"""
function evaluate!(dihedralsCos::Array{Forcefield.DihedralCos}, state::Common.State;
    do_forces = false)

    v12 = zeros(Float64, 3)
    v32 = zeros(Float64, 3)
    v34 = zeros(Float64, 3)
    f1 = zeros(Float64, 3)
    f3 = zeros(Float64, 3)
    f4 = zeros(Float64, 3)
    m = zeros(Float64, 3)
    n = zeros(Float64, 3)
    
    energy = 0.0

    for dihedral in dihedralsCos
        
        @views @. v12 = state.xyz[dihedral.a2, :] - state.xyz[dihedral.a1, :]
        @views @. v32 = state.xyz[dihedral.a2, :] - state.xyz[dihedral.a3, :]
        @views @. v34 = state.xyz[dihedral.a4, :] - state.xyz[dihedral.a3, :]
        m = cross(v12, v32)
        n = cross(v32, v34)
        d32Sq = dot(v32,v32)
        d32 = sqrt(d32Sq)
        phi = atan(d32 * dot(v12, n), dot(m, n))# - pi
        
        energy += dihedral.k * (1.0 + cos(dihedral.mult * phi - dihedral.θ))
        
        if do_forces
            dVdphi_x_d32 = dihedral.k * dihedral.mult * sin(dihedral.θ - dihedral.mult * phi) * d32
            @. f1 = m * (-dVdphi_x_d32 / dot(m, m))
            @. f4 = n * ( dVdphi_x_d32 / dot(n, n))
            # @. f3 = f4 * (dot(v34, v32)/d32Sq - 1.0) - f1 * (dot(v12, v32)/d32Sq)
            @. f3 = -f4
            @. f3 -= f1 * (dot(v12, v32)/d32Sq)
            @. f3 += f4 * (dot(v34, v32)/d32Sq)
            
            @. state.forces[dihedral.a1, :] -= f1
            @. state.forces[dihedral.a2, :] -= (-f1 - f3 - f4)
            @. state.forces[dihedral.a3, :] -= f3
            @. state.forces[dihedral.a4, :] -= f4
        end
    end

    state.energy.eDihedral = energy

    return energy
end

@doc raw"""
    evaluate!(atoms::Array{Forcefield.Atom}, state::Common.State, do_forces::Bool = false)::Float64

Evaluate an array of [Forcefield.Components](#Components-1) using the current [`Common.State`](@ref),
calculate and update state.energy according to the equations defined in each component.
If `do_forces` flag is set to `true`, calculate and update `state.forces`.
Non-bonded interactions are only assessed if the distance between atoms is below the defined `cut_off` value.
Return the component energy value (kJ mol⁻¹).

# Examples
```julia-repl
julia> Forcefield.evaluate!(bonds, state)
0.500
```

See also: [`evalenergy!`](@ref) [`Forcefield.HarmonicBond`](@ref) [`Forcefield.HarmonicAngle`](@ref)
[`Forcefield.DihedralCos`](@ref) [`Forcefield.Atom`](@ref)
"""
function evaluate!(atoms::Array{Forcefield.Atom}, state::Common.State;
    do_forces::Bool = false, cut_off::Float64 = 2.0)

    eLJ = 0.0
    eLJ14 = 0.0
    eCoulomb = 0.0
    eCoulomb14 = 0.0

    n_atoms = length(atoms)
    vij = zeros(Float64, 3)
    cut_offSq = cut_off*cut_off
    exclude_idx::Int64 = 1
    exclude::Int64 = 1
    
    #Calculate nonbonded interactions
    for i in 1:(n_atoms - 1)
        atomi::Forcefield.Atom = atoms[i]
        
        # set the exclution index to the correct location and extract the exclude atom index
        exclude_idx = 1
        while atomi.excls[exclude_idx] <= i && exclude_idx < length(atomi.excls)
            exclude_idx += 1
        end
        exclude = atomi.excls[exclude_idx]

        for j in (i+1):(n_atoms)
            
            if j == exclude
                exclude_idx += 1
                exclude = atomi.excls[exclude_idx]
                continue
            end
            atomj::Forcefield.Atom = atoms[j]
            
            @views @. vij = state.xyz[j, :] - state.xyz[i, :]

            #Check if the distance between the two atoms is below cut-off
            dijSq = dot(vij, vij)
            if dijSq > cut_offSq
                continue
            end

            #Calculate energy (σ and ϵ already have the necessary constants multiplied)
            sij = atomi.σ + atomj.σ
            eij = atomi.ϵ * atomj.ϵ
            lj6 = (sij * sij/dijSq) ^ 3
            eLJ += eij * (lj6 * lj6 - lj6)
            ecoul = atomi.q * atomj.q / sqrt(dijSq)
            eCoulomb += ecoul

            #Calculate forces, if requested
            if do_forces
                fc = (24.0 * eij * (lj6 - 2.0 * lj6 * lj6) - ecoul) / dijSq
                @. vij *= fc
                @. state.forces[i, :] += vij
                @. state.forces[j, :] -= vij
            end
        end
    end

    eLJ *= 4.0
    state.energy.eLJ = eLJ
    state.energy.eCoulomb = eCoulomb

    #Calculate 1-4 interactions
    evdw_scale = 0.5
    ecoul_scale = 0.833333

    for i in 1:(n_atoms)
        atomi = atoms[i]
        for j in atomi.pairs
            
            atomj = atoms[j]
            @views @. vij = state.xyz[j, :] - state.xyz[i, :]
            dijSq = dot(vij, vij)
            
            sij = atomi.σ + atomj.σ
            eij = evdw_scale * atomi.ϵ * atomj.ϵ
            qij = ecoul_scale * atomi.q * atomj.q
            lj6 = (sij * sij/dijSq) ^ 3
            eLJ14 += eij * (lj6 * lj6 - lj6)
            ecoul = qij / sqrt(dijSq)
            eCoulomb14 += ecoul

            # Calculate forces, if requested
            if do_forces
                fc = (24.0 * eij * (lj6 - 2.0 * lj6 * lj6) - ecoul) / dijSq
                @. vij *= fc
                @. state.forces[i, :] += vij
                @. state.forces[j, :] -= vij
            end
        end
    end

    eLJ14 *= 4.0
    state.energy.eLJ14 = eLJ14
    state.energy.eCoulomb14 = eCoulomb14

    return eLJ + eLJ14 + eCoulomb + eCoulomb14
end

@doc raw"""
    evalenergy!(topology::Forcefield.Topology, state::Common.State[, cut_off::Float64 = 2.0, do_forces::Bool = false])::Float64

Evaluate the current [`Common.State`](@ref) energy according to the defined [`Forcefield.Topology`](@ref).
If `do_forces` bool is set to `true`, calculate and update `state.forces`.
Non-bonded interactions are only assessed if the distance between atoms is below the defined `cut_off` value.
Return `state.energy.eTotal` value (kJ mol⁻¹).

# Examples
```julia-repl
julia> Forcefield.evalenergy!(topology, state, cut_off = Inf)
0.500
```

See also: [`evaluate!`](@ref)
"""
function evalenergy!(topology::Forcefield.Topology, state::Common.State;
    cut_off::Float64 = 2.0, do_forces = false)

    energy = 0.0
    energy += evaluate!(topology.bonds, state, do_forces = do_forces)
    energy += evaluate!(topology.angles, state, do_forces = do_forces)
    energy += evaluate!(topology.atoms, state, do_forces = do_forces, cut_off = cut_off)
    # energy += evaldihedralRB(topology.dihedralsRB, state, do_forces = do_forces)
    energy += evaluate!(topology.dihedralsCos, state, do_forces = do_forces)
    return energy
end