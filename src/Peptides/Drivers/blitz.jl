using Random
using ProtoSyn.Calculators: EnergyFunction
using ProtoSyn: AbstractSelection, getdihedral, promote
using ProtoSyn.Peptides: Rotamers, Dihedral
using ProtoSyn.Drivers: Driver, DriverState

Base.@kwdef mutable struct RotamerBlitzState <: DriverState
    step::Int        = 0
    converged::Bool  = false
    completed::Bool  = false
    stalled::Bool    = false
    acceptance_count = 0
end


mutable struct RotamerBlitz <: Driver
    eval!::Union{Function, EnergyFunction}
    rotamer_library::Dict{String, ProtoSyn.Peptides.Rotamers.BBD_RotamerLibrary}
    n_first::Int
    max_steps::Int
    selection::Opt{AbstractSelection}
end

RotamerBlitz(eval!::Union{Function, EnergyFunction}, rotamer_library::Dict{String, ProtoSyn.Peptides.Rotamers.BBD_RotamerLibrary}, n_first::Int, max_steps::Int) = begin
    return RotamerBlitz(eval!, rotamer_library, n_first, max_steps, nothing)
end

function (driver::RotamerBlitz)(cb::Opt{F}, pose::Pose) where {F<:Function}
    
    T = eltype(pose.state)
    driver_state = RotamerBlitzState()
    if driver.selection === nothing
        residues = collect(eachresidue(pose.graph))
    else
        residues = promote(driver.selection, Residue)(pose, gather = true)
    end

    cb !== nothing && cb(pose, driver_state)

    for step in 1:driver.max_steps
        for residue in shuffle(residues)
            
            # Residue results hold the energy of each rotamer conformation
            # The best rotamer conformation will be applied to this residue
            energy          = driver.eval!(pose)
            rotamer         = Rotamers.get_rotamer(pose, residue)
            residue_results = Dict{Rotamers.Rotamer, T}(rotamer => energy)

            # Get the correct rotamer stack
            phi = getdihedral(pose.state, Dihedral.phi(residue))
            psi = getdihedral(pose.state, Dihedral.phi(residue))
            !(residue.name in keys(driver.rotamer_library)) && continue
            rotamer_stack = driver.rotamer_library[residue.name][phi, psi]

            # Iterate the 'n_first' rotamers, saving the evaluated energy
            for rotamer in rotamer_stack.rotamers[1:driver.n_first]
                # Note: Rotamers.apply! doesn't sync! the pose
                exact_rotamer = Rotamers.apply!(pose.state, rotamer, residue)
                sync!(pose)
                energy = driver.eval!(pose)
                residue_results[exact_rotamer] = energy
            end

            # Re-apply the best conformation (exactly)
            min_energy, rotamer = findmin(residue_results)
            Rotamers.apply!(pose.state, rotamer, residue)
        end

        driver_state.step += 1
        cb !== nothing && cb(pose, driver_state) 
    end
    
    driver_state.completed = true
    driver_state
end

(driver::RotamerBlitz)(pose::Pose) = driver(nothing, pose)