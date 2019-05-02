using ProtoSyn

@doc raw"""
    evaluate!(bonds::Vector{Forcefield.Amber.HarmonicBond}, state::Common.State[, do_forces::Bool = false])::Float64
"""
function evaluate!(bonds::Vector{Forcefield.Amber.HarmonicBond}, state::Common.State; do_forces::Bool = false)::Float64

    energy = 0.0
    Δ      = 0.0
    d12    = 0.0
    v12    = [0.0, 0.0, 0.0]
    forces = state.forces
    coords = state.xyz

    for bond in bonds
        d12 = 0.0
        @inbounds for i=1:3
            Δ = coords[bond.a2, i] - coords[bond.a1, i]
            v12[i] = Δ
            d12 += Δ ^ 2
        end
        
        d12 = sqrt(d12)
        dr = d12 - bond.b0
        energy += bond.k * dr * dr

        if do_forces
            Δ = bond.k * dr / d12
            @inbounds for i=1:3
                forces[bond.a1, i] += Δ * v12[i]
                forces[bond.a2, i] -= Δ * v12[i]
            end
        end
    end

    energy *= .5
    Common.set_energy_component(state.energy, :bond, energy)
    return energy
end

@doc raw"""
    evaluate!(angles::Vector{Forcefield.Amber.HarmonicAngle}, state::Common.State, do_forces::Bool = false)::Float64
"""
function evaluate!(angles::Vector{HarmonicAngle}, state::Common.State; do_forces::Bool = false)::Float64

    energy     = 0.0
    delta      = 0.0
    d12        = 0.0
    v12        = [0.0, 0.0, 0.0]
    v32        = [0.0, 0.0, 0.0]
    forces     = state.forces
    coords     = state.xyz

    @inbounds for angle in angles
        ctheta = 0.0
        d12Sq  = 0.0
        d32Sq  = 0.0

        @inbounds for i=1:3
            a1 = coords[angle.a1, i]
            a2 = coords[angle.a2, i]
            a3 = coords[angle.a3, i]

            delta  = a2 - a1
            v12[i] = delta
            d12Sq += delta*delta

            delta2  = a2 - a3
            v32[i] = delta2
            d32Sq += delta2*delta2

            ctheta += delta*delta2
        end

        d12xd32 = sqrt(d12Sq * d32Sq)
        ctheta /= d12xd32
        dtheta  = acos(ctheta) - angle.θ
        kdtheta = angle.k * dtheta
        energy += kdtheta * dtheta
        
        if do_forces
            fc = kdtheta / sqrt(1.0 - ctheta * ctheta)
            @inbounds for i=1:3
                f1 = fc*(ctheta * v12[i]/d12Sq - v32[i]/d12xd32)
                f3 = fc*(ctheta * v32[i]/d32Sq - v12[i]/d12xd32)

                forces[angle.a1, i] += f1 
                forces[angle.a3, i] += f3
                forces[angle.a2, i] -= (f3+f1)
            end
        end
    end

    energy *= 0.5
    state.energy.comp["eAngle"] = energy
    state.energy.eTotal = energy
    return energy
end

@doc raw"""
    evaluate!(dihedralsCos::Vector{Forcefield.Amber.DihedralCos}, state::Common.State, do_forces::Bool = false)::Float64
"""
function evaluate!(dihedralsCos::Vector{DihedralCos}, state::Common.State; do_forces = false)::Float64

    energy  = 0.0
    d12 = 0.0
    d32 = 0.0
    d34 = 0.0
    v12 = [0.0, 0.0, 0.0]
    v32 = [0.0, 0.0, 0.0]
    v34 = [0.0, 0.0, 0.0]
    m   = [0.0, 0.0, 0.0]
    n   = [0.0, 0.0, 0.0]
    forces = state.forces
    coords = state.xyz

    for dihedral in dihedralsCos
        d32Sq  = 0.0
        d3432  = 0.0
        d1232  = 0.0

        @inbounds for i=1:3
            d12 = coords[dihedral.a2, i] - coords[dihedral.a1, i]
            d32 = coords[dihedral.a2, i] - coords[dihedral.a3, i]
            d34 = coords[dihedral.a4, i] - coords[dihedral.a3, i]
            v12[i]  = d12

            v32[i]  = d32
            d32Sq  += d32 ^ 2

            v34[i]  = d34

            d3432  += d34 * d32
            d1232  += d12 * d32
        end
        
        # m = v12 × v32 = cross(v12, v32)
        m[1] = v12[2]*v32[3] - v12[3]*v32[2]
        m[2] = v12[3]*v32[1] - v12[1]*v32[3]
        m[3] = v12[1]*v32[2] - v12[2]*v32[1]
        
        # n = v32 × v34 = cross(v32, v34)  
        n[1] = v32[2]*v34[3] - v32[3]*v34[2]
        n[2] = v32[3]*v34[1] - v32[1]*v34[3]
        n[3] = v32[1]*v34[2] - v32[2]*v34[1]

        d32     = sqrt(d32Sq)
        phi     = atan(d32 * dot(v12, n), dot(m, n))
        energy += dihedral.k * (1.0 + cos(dihedral.mult * phi - dihedral.θ))
        
        if do_forces
            dVdphi_x_d32 = dihedral.k * dihedral.mult * sin(dihedral.θ - dihedral.mult * phi) * d32
            # Constants:
            c1   = -dVdphi_x_d32 / dot(m, m)
            c4   = dVdphi_x_d32 / dot(n, n)
            c3_1 = (d3432/d32Sq - 1.0)
            c3_2 = (d1232/d32Sq)
            @inbounds for i=1:3
                f1 = m[i] * c1
                f4 = n[i] * c4
                f3 = f4 * c3_1 - f1 * c3_2

                forces[dihedral.a1, i] -= f1
                forces[dihedral.a2, i] -= (-f1 - f3 - f4)
                forces[dihedral.a3, i] -= f3
                forces[dihedral.a4, i] -= f4
            end
        end
    end

    state.energy.comp["eDihedral"] = energy
    state.energy.eTotal = energy
    return energy
end

@doc raw"""
    evaluate!(atoms::Vector{Forcefield.Amber.Atom}, state::Common.State[, do_forces::Bool = false, cut_off::Float64 = 2.0])::Float64

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
function evaluate!(atoms::Vector{Atom}, state::Common.State; do_forces::Bool = false, cut_off::Float64 = 2.0)::Float64

    eLJ         = 0.0
    eLJ14       = 0.0
    eCoulomb    = 0.0
    eCoulomb14  = 0.0
    vij         = [0.0, 0.0, 0.0]
    exclude_idx = 1
    exclude     = 1
    cut_offSq   = cut_off ^ 2
    coords      = state.xyz
    forces      = state.forces
    n_atoms     = length(atoms)
    
    #Calculate nonbonded interactions
    
    for i in 1:(n_atoms - 1)
        @inbounds atomi = atoms[i]
        
        # set the exclution index to the correct location and extract the exclude atom index
        len_exclusions = length(atomi.excls)
        if length(atomi.excls) > 0
            exclude_idx = 1
            @inbounds while atomi.excls[exclude_idx] <= i && exclude_idx < len_exclusions
                exclude_idx += 1
            end
            @inbounds exclude = atomi.excls[exclude_idx]
        end
        
        # for j in (i+1):(n_atoms)
        #     if j == exclude
        #         exclude_idx += 1
        #         @inbounds exclude = atomi.excls[exclude_idx]
        #         continue
        #     end
        
        ptr = state.nb.pointer[i]
        while state.nb.list[ptr] > 0
            @inbounds j = state.nb.list[ptr]
            ptr += 1
            
            if j == exclude
                exclude_idx += 1
                @inbounds exclude = atomi.excls[exclude_idx]
                continue
            end


            @inbounds atomj = atoms[j]
            
            dijSq = 0.0
            @inbounds for k=1:3
                deltaij = coords[j, k] - coords[i, k]
                dijSq  += deltaij * deltaij
                vij[k]  = deltaij
            end

            #Check if the distance between the two atoms is below cut-off
            if dijSq > cut_offSq
                continue
            end

            #Calculate energy (σ and ϵ already have the necessary constants multiplied)
            sij = atomi.σ + atomj.σ
            eij = atomi.ϵ * atomj.ϵ
            lj6 = (sij * sij/dijSq) ^ 3
            eLJ += eij * (lj6 * lj6 - lj6)
            # eCoulomb = eCoulomb_λ * atomi.q * atomj.q / sqrt(dijSq)
            eCoulomb = atomi.q * atomj.q / sqrt(dijSq)

            #Calculate forces, if requested
            if do_forces
                fc = (24.0 * eij * (lj6 - 2.0 * lj6 * lj6) - eCoulomb) / dijSq
                @inbounds for k=1:3
                    t = vij[k] * fc
                    forces[i, k] += t
                    forces[j, k] -= t
                end
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
            dijSq = 0.0
            @inbounds for k=1:3
                deltaij = coords[j, k] - coords[i, k]
                vij[k]  = deltaij
                dijSq  += deltaij * deltaij
            end
            
            sij = atomi.σ + atomj.σ
            eij = evdw_scale * atomi.ϵ * atomj.ϵ 
            qij = ecoul_scale * atomi.q * atomj.q
            lj6 = (sij * sij/dijSq) ^ 3
            eLJ14 += eij * ((lj6 ^ 2) - lj6)
            # eCoulomb14 = eCoulomb_λ * qij / sqrt(dijSq)
            eCoulomb14 = qij / sqrt(dijSq)

            # Calculate forces, if requested
            if do_forces
                fc = (24.0 * eij * (lj6 - 2.0 * lj6 * lj6) - eCoulomb14) / dijSq
                @inbounds for k=1:3
                    t = vij[k] * fc
                    forces[i, k] += t
                    forces[j, k] -= t
                end
            end
        end
    end
    
    eLJ14 *= 4.0
    state.energy.comp["eLJ14"] = eLJ14
    state.energy.comp["eCoulomb14"] = eCoulomb14

    energy = eLJ + eLJ14 + eCoulomb + eCoulomb14
    state.energy.eTotal = energy
    return energy
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
function evaluate!(topology::Topology, state::Common.State; cut_off::Float64 = 1.2, do_forces::Bool = false)::Float64
    
    energy =  evaluate!(topology.bonds, state, do_forces = do_forces)
    energy += evaluate!(topology.angles, state, do_forces = do_forces)
    energy += evaluate_OLD!(topology.atoms, state, do_forces = do_forces, cut_off = cut_off)
    energy += evaluate!(topology.dihedralsCos, state, do_forces = do_forces)
    state.energy.comp["Amber"] = energy
    state.energy.eTotal = energy
    return energy
end

#! DEPRECATED
function evaluate_OLD!(atoms::Vector{Atom}, state::Common.State; do_forces::Bool = false, cut_off::Float64 = 2.0,
    eCoulomb_λ::Float64 = 1.0)::Float64

    eLJ         = .0
    eLJ14       = .0
    eCoulomb    = .0
    eCoulomb14  = .0
    n_atoms     = length(atoms)
    vij         = [.0, .0, .0]
    cut_offSq   = cut_off*cut_off
    exclude_idx = 1
    exclude     = 1
    coords      = state.xyz
    forces      = state.forces
    
    #Calculate nonbonded interactions
    for i in 1:(n_atoms - 1)
        atomi = atoms[i]
        
        # set the exclution index to the correct location and extract the exclude atom index
        if length(atomi.excls) > 0
            exclude_idx = 1
            while atomi.excls[exclude_idx] <= i && exclude_idx < length(atomi.excls)
                exclude_idx += 1
            end
            exclude = atomi.excls[exclude_idx]
        end

        for j in (i+1):(n_atoms)
        
            if j == exclude
                exclude_idx += 1
                exclude = atomi.excls[exclude_idx]
                continue
            end
            atomj = atoms[j]
            
            dijSq = .0
            @inbounds for k=1:3
                deltaij = coords[j, k] - coords[i, k]
                vij[k]  = deltaij
                dijSq  += deltaij * deltaij
            end

            #Check if the distance between the two atoms is below cut-off
            if dijSq > cut_offSq
                continue
            end

            #Calculate energy (σ and ϵ already have the necessary constants multiplied)
            sij = atomi.σ + atomj.σ
            eij = atomi.ϵ * atomj.ϵ
            lj6 = (sij * sij/dijSq) ^ 3
            eLJ += eij * (lj6 * lj6 - lj6)
            ecoul = atomi.q * atomj.q / sqrt(dijSq)
            eCoulomb += eCoulomb_λ * ecoul

            #Calculate forces, if requested
            if do_forces
                fc = (24.0 * eij * (lj6 - 2.0 * lj6 * lj6) - ecoul) / dijSq
                @inbounds for k=1:3
                    t = vij[k] * fc
                    forces[i, k] += t
                    forces[j, k] -= t
                end
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
            dijSq = .0
            @inbounds for k=1:3
                deltaij = coords[j, k] - coords[i, k]
                vij[k]  = deltaij
                dijSq  += deltaij * deltaij
            end
            
            sij = atomi.σ + atomj.σ
            eij = evdw_scale * atomi.ϵ * atomj.ϵ 
            qij = ecoul_scale * atomi.q * atomj.q
            lj6 = (sij * sij/dijSq) ^ 3
            eLJ14 += eij * ((lj6 ^ 2) - lj6)
            ecoul = qij / sqrt(dijSq)
            eCoulomb14 += eCoulomb_λ * ecoul

            # Calculate forces, if requested
            if do_forces
                fc = (24.0 * eij * (lj6 - 2.0 * lj6 * lj6) - ecoul) / dijSq
                @inbounds for k=1:3
                    t = vij[k] * fc
                    forces[i, k] += t
                    forces[j, k] -= t
                end
            end
        end
    end

    eLJ14 *= 4.0
    state.energy.comp["eLJ14"] = eLJ14
    state.energy.comp["eCoulomb14"] = eCoulomb14

    energy = eLJ + eLJ14 + eCoulomb + eCoulomb14
    state.energy.eTotal = energy
    return energy
end