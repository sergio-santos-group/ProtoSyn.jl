using ProtoSyn
using Printf
using LinearAlgebra

include("ilsrr_config.jl")

print_status = @Common.callback 1000 function cb_status(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.Driver, max_force::Float64, gamma::Float64, args...)
    Print.status(@sprintf("(%5s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Max Force: %10.3e | Gamma: %10.3e\n", "SD", step, st.energy.eTotal, st.energy.comp["amber"],
        st.energy.comp["other"], max_force, gamma), stdout)
end

print_structure = @Common.callback 1000 function cb_structure(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.Driver, args...)
    Print.as_pdb(xyz_destination, state, metadata)
end

function amber_evaluator!(st::Common.State, do_forces::Bool)
    fill!(st.forces, zero(Float64))
    e  = Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
    # e += Forcefield.Restraints.evaluate!(dihedral_restraints, st, do_forces=do_forces)
    st.energy.eTotal = e
    return e
end

minimizer   = Drivers.SteepestDescent.Driver(amber_evaluator!, 1000, 100.0, 1e-4, print_status)
loop_closer = Drivers.SteepestDescent.Driver(amber_evaluator!, 10000, 100.0, 1e-4, print_status, print_structure)

state, metadata     = Common.load_from_pdb("saved/4500.pdb")
metadata.blocks     = Common.compile_blocks(metadata.residues, ss)
metadata.ss         = Common.compile_ss(metadata.residues, ss)
amber_topology      = Forcefield.Amber.load_from_json(input_amber_json)
dihedral_restraints = Forcefield.Restraints.lock_block_bb(metadata, fbw = dihedral_fb_width, k = 1e10)
blockrot_mutator    = Mutators.Blockrot.BlockrotMutator(metadata.blocks, () -> (randn() * blockrot_mutator.step_size), 1.0, π/4, loop_closer)

# println("Starting initial minimization")
# loop_closer.run!(state, minimizer)
# println("Finished initial minimization")
Print.as_pdb(xyz_destination, state, metadata)
r_count = Mutators.Blockrot.run!(state, blockrot_mutator)
Print.as_pdb(xyz_destination, state, metadata)