mutable struct PlainCutoff <: IForcefieldConfig
    cutoff::Float64
    PlainCutoff(c::Float64) = begin
        if c < 0.0
            error("cutoff must be a non-negative value")
        end
        new(c)
    end
end



# proxy for dispatching the correct function
# Calculators.eval!(s::State, c::AtomContainer, f::Bool=true) = 
Calculators.eval!(s::State, c::ComponentContainer{Atom}, f::Bool=true) = 
    Calculators.eval!(s, c.items, f, c.config)


function Calculators.eval!(state::State, atoms::Vector{Atom}, do_forces::Bool, ::Nothing)
    
    eLJ = 0.0
    eCoulomb = 0.0
    forces = state.forces
    coords = state.coords
    natoms = state.size

    @inbounds for i = 1:natoms
        atom_i = atoms[i]
        type_i = atom_i.type
        
        exclude = atom_i.exclusions===nothing ? -1 : atom_i.exclusions[1]
        excl_idx = 2

        for j = (i+1):natoms
            
            if j == exclude
                exclude = atom_i.exclusions[excl_idx]
                excl_idx += 1
                continue
            end

            atom_j = atoms[j]
            type_j = atom_j.type

            @nexprs 3 u -> v12_u = coords[j,u] - coords[i,u]

            d12sq = @dot u v12_u v12_u

            # lennard-jones potential
            σ12 = type_i.σ + type_j.σ
            ϵ12 = type_i.ϵ * type_j.ϵ

            lj2 = σ12*σ12/d12sq
            lj6 = lj2*lj2*lj2
            eLJ += ϵ12 * (lj6 * lj6 - lj6)

            # coulomb potential
            ecoul = 138.935458 * atom_i.q * atom_j.q / sqrt(d12sq)
            eCoulomb += ecoul

            if do_forces
                fc = (24.0 * ϵ12 * (lj6 - 2.0 * lj6 * lj6) - ecoul) / d12sq
                @nexprs 3 u -> begin
                    forces[i, u] += fc * v12_u
                    forces[j, u] -= fc * v12_u
                end
            end
        end

    end
    
    state.energy.LJ = 4.0 * eLJ
    state.energy.Coulomb = eCoulomb
    eCoulomb + 4.0 * eLJ
end



function Calculators.eval!(state::State, atoms::Vector{Atom}, do_forces::Bool, cfg::PlainCutoff)
    # println("PlainCutoff")

    eLJ = 0.0
    eCoulomb = 0.0
    forces = state.forces
    coords = state.coords
    natoms = state.size

    cutoff_sq = cfg.cutoff^2
    #nattemps = 0
    #naccepts = 0
    @inbounds for i = 1:natoms
        atom_i = atoms[i]
        type_i = atom_i.type
        
        exclude = atom_i.exclusions===nothing ? -1 : atom_i.exclusions[1]
        excl_idx = 2

        @nexprs 3 u -> fi_u = 0.0

        ptr = state.pairlist.offset[i]
        while state.pairlist.list[ptr] > 0
           j = state.pairlist.list[ptr]
           ptr += 1
        # for j = (i+1):natoms
            
            if j == exclude
                exclude = atom_i.exclusions[excl_idx]
                excl_idx += 1
                continue
            end

            atom_j = atoms[j]
            type_j = atom_j.type

            @nexprs 3 u -> v12_u = coords[j,u] - coords[i,u]

            d12sq = @dot u v12_u v12_u
            #nattemps += 1
            (d12sq > cutoff_sq) && continue
            #naccepts += 1
            
            # lennard-jones potential
            σ12 = type_i.σ + type_j.σ
            ϵ12 = type_i.ϵ * type_j.ϵ

            lj2 = σ12*σ12/d12sq
            lj6 = lj2*lj2*lj2
            eLJ += ϵ12 * (lj6 * lj6 - lj6)

            #elj = ϵ12 * (lj6 * lj6 - lj6)
            #if elj > 0.0
            #    println("  $i - $j: $elj")
            #end
            
            
            # coulomb potential
            ecoul = 138.935458 * atom_i.q * atom_j.q / sqrt(d12sq)
            eCoulomb += ecoul

            if do_forces
                fc = (24.0 * ϵ12 * (lj6 - 2.0 * lj6 * lj6) - ecoul) / d12sq
                @nexprs 3 u -> begin
                    # forces[i, u] += fc * v12_u
                    fi_u += fc * v12_u
                    forces[j, u] -= fc * v12_u
                end
            end
        end
        @nexprs 3 u -> forces[i, u] += fi_u

        # println("move update force[i] here")

    end
    #println("nattemps = $nattemps naccepts = $naccepts")
    state.energy.LJ = 4.0 * eLJ
    state.energy.Coulomb = eCoulomb
    eCoulomb + 4.0 * eLJ
end