@doc raw"""
    evaluate!(bonds::Array{Forcefield.HarmonicBond}, state::Common.State[, do_forces::Bool = false])::Float64
"""
function evaluate!(bonds::Vector{HarmonicBond}, state::Common.State;
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
    state.energy.comp["eBond"] = 0.5 * energy
    state.energy.eTotal = 0.5 * energy
    return 0.5 * energy
end

@doc raw"""
    evaluate!(angles::Array{Forcefield.HarmonicAngle}, state::Common.State, do_forces::Bool = false)::Float64
"""
function evaluate!(angles::Vector{HarmonicAngle}, state::Common.State;
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
    state.energy.comp["eAngle"] = 0.5 * energy
    state.energy.eTotal = 0.5 * energy
    return 0.5 * energy
end

@doc raw"""
    evaluate!(dihedralsCos::Array{Forcefield.DihedralCos}, state::Common.State, do_forces::Bool = false)::Float64
"""
function evaluate!(dihedralsCos::Vector{DihedralCos}, state::Common.State;
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
        phi = atan(d32 * dot(v12, n), dot(m, n))
        println(dihedral, " ", rad2deg(phi))
        
        energy += dihedral.k * (1.0 + cos(dihedral.mult * phi - dihedral.θ))
        
        if do_forces
            dVdphi_x_d32 = dihedral.k * dihedral.mult * sin(dihedral.θ - dihedral.mult * phi) * d32
            f1 .= m .* (-dVdphi_x_d32 / dot(m, m))
            f4 .= n .* ( dVdphi_x_d32 / dot(n, n))
            f3 .= f4 .* (dot(v34, v32)/d32Sq - 1.0) .- f1 .* (dot(v12, v32)/d32Sq)
            
            # @. f3 = -f4
            # @. f3 -= f1 * (dot(v12, v32)/d32Sq)
            # @. f3 += f4 * (dot(v34, v32)/d32Sq)
            
            @. state.forces[dihedral.a1, :] -= f1
            @. state.forces[dihedral.a2, :] -= (-f1 - f3 - f4)
            @. state.forces[dihedral.a3, :] -= f3
            @. state.forces[dihedral.a4, :] -= f4
        end
    end
    state.energy.comp["eDihedral"] = energy
    state.energy.eTotal = energy
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
julia> Forcefield.Amber.evaluate!(bonds, state)
0.500
```

See also: [`evaluate!`](@ref) [`Amber.HarmonicBond`](@ref Forcefield) [`Amber.HarmonicAngle`](@ref Forcefield)
[`Amber.DihedralCos`](@ref Forcefield) [`Amber.Atom`](@ref Forcefield)
"""
function evaluate!(atoms::Vector{Atom}, state::Common.State;
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
        atomi::Atom = atoms[i]
        
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
            atomj::Atom = atoms[j]
            
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
    state.energy.comp["eLJ"] = eLJ
    state.energy.comp["eCoulomb"] = eCoulomb

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
    state.energy.comp["eLJ14"] = eLJ14
    state.energy.comp["eCoulomb14"] = eCoulomb14

    state.energy.eTotal = (eLJ + eLJ14 + eCoulomb + eCoulomb14)
    return eLJ + eLJ14 + eCoulomb + eCoulomb14
end

@doc raw"""
    evaluate!(topology::Forcefield.Topology, state::Common.State[, cut_off::Float64 = 2.0, do_forces::Bool = false])::Float64

Evaluate the current [`Common.State`](@ref) energy according to the defined [`Amber.Topology`](@ref Forcefield).
If `do_forces` bool is set to `true`, calculate and update `state.forces`.
Non-bonded interactions are only assessed if the distance between atoms is below the defined `cut_off` value.
Return `state.energy.eTotal` value (kJ mol⁻¹).

# Examples
```julia-repl
julia> Forcefield.Amber.evaluate!(topology, state, cut_off = Inf)
0.500
```

See also: [`Amber.evaluate!`](@ref)
"""
function evaluate!(topology::Topology, state::Common.State;
    cut_off::Float64 = 2.0, do_forces = false)
    
    energy =  evaluate!(topology.bonds, state, do_forces = do_forces)
    energy += evaluate!(topology.angles, state, do_forces = do_forces)
    energy += evaluate!(topology.atoms, state, do_forces = do_forces, cut_off = cut_off)
    energy += evaluate!(topology.dihedralsCos, state, do_forces = do_forces)
    state.energy.eTotal = energy
    return energy
end