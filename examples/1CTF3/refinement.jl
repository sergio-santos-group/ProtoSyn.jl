using ProtoSyn
using Printf

state, metadata     = Common.load_from_pdb("data/1ctf_native.pdb")
amber_topology      = Forcefield.Amber.load_from_json("data/1ctf_amber_top.json")
metadata.ss         = Common.compile_ss(metadata.dihedrals, "CEEEEEEECCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCEEEEEEECHHHHHHHHHHHHHHCCEEEEC")
metadata.blocks     = Common.compile_blocks(metadata.residues, "CEEEEEEECCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCEEEEEEECHHHHHHHHHHHHHHCCEEEEC")
contact_restraints  = Forcefield.Restraints.load_distance_restraints_from_file("data/contact_map_raptorx.txt", metadata, k = 10.0, threshold = 0.31)
dihedral_restraints = Forcefield.Restraints.lock_block_bb(metadata, fbw = 10.0, k = 100.0)
nb_dhs              = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dhs          = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dhs)

const out_pdb       = open("out/refinement.pdb", "w")
const out_pdb_best  = open("out/refinement_best.pdb", "w")
const out_ene       = open("out/refinement.dat", "w")

print_status_init = @Common.callback 10 function _print_status_init(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.Driver, max_force::Float64, gamma::Float64, args...)
    Print.status(@sprintf("(%5s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Max Force: %10.3e | Gamma: %10.3e\n", "I_SD", step, st.energy.eTotal,
        st.energy.comp["amber"], st.energy.comp["other"], max_force, gamma), color = :red)
end

print_status_refnm = @Common.callback 1 function _print_status_refnm(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.Driver, ar::Float64, args...)
    Print.status(@sprintf("(%5s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | DH: %8.2e | CS: %8.2e | BR: %8.2f | AR: %.4f | T: %4.2f\n",
        "REFNM", step, st.energy.eTotal, st.energy.comp["amber"], st.energy.comp["other"], soft_dh_mutator.step_size, soft_cs_mutator.step_size, soft_br_mutator.step_size,
        ar, dr.temperature), color = :green)
end

print_status_loop = @Common.callback 100 function _print_status_loop(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.Driver, max_force::Float64, gamma::Float64, args...)
    Print.status(@sprintf("⤷(%4s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Max Force: %10.3e | Gamma: %10.3e\n", "LOOP", step, st.energy.eTotal,
        st.energy.comp["amber"], st.energy.comp["other"], max_force, gamma))
end

print_status_smpl = @Common.callback 100 function _print_status_smpl(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.Driver, max_force::Float64, gamma::Float64, args...)
    Print.status(@sprintf("⤷(%4s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Max Force: %10.3e | Gamma: %10.3e\n", "SMPL", step, st.energy.eTotal,
        st.energy.comp["amber"], st.energy.comp["other"], max_force, gamma))
end

print_structure_rfnm = @Common.callback 1 function _print_structure_rfnm(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    Print.as_pdb(out_pdb, st, metadata, step = step)
    flush(out_pdb)
end

outer_best_energy = Inf
print_structure_best = @Common.callback 1 function _print_structure_best(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    global outer_best_energy
    if st.energy.eTotal < outer_best_energy
        Print.as_pdb(out_pdb_best, st, metadata, step = step)
        flush(out_pdb_best)
        printstyled(@sprintf("(ILSRR) New outer best defined: ⚡E: %10.3e (old) ▶️ %10.3e (new)\n", outer_best_energy, st.energy.eTotal), color = :red)
        outer_best_energy = st.energy.eTotal
    end
end

function amber_cr_dr_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    amber_e    = Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
    contact_e  = Forcefield.Restraints.evaluate!(contact_restraints, st, do_forces=do_forces)
    dihedral_e = Forcefield.Restraints.evaluate!(dihedral_restraints, st, do_forces=do_forces)
    energy = amber_e + contact_e + dihedral_e
    st.energy.comp["other"] = contact_e + dihedral_e
    st.energy.eTotal = energy
    return energy
end

function amber_cr_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    amber_e   = Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
    contact_e = Forcefield.Restraints.evaluate!(contact_restraints, st, do_forces=do_forces)
    energy    = amber_e + contact_e
    st.energy.comp["other"] = contact_e
    st.energy.eTotal = energy
    return energy
end

function amber_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    return Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
end

initial_minimizer = Drivers.SteepestDescent.Driver(amber_evaluator!, 10000, 2000.0, 0.01, true, print_status_init)
sampler_minimizer = Drivers.SteepestDescent.Driver(amber_cr_dr_evaluator!, 10000, 2000.0, 0.01, true, print_status_smpl)
loop_closer       = Drivers.SteepestDescent.Driver(amber_evaluator!, 10000, 2000.0, 0.01, true, print_status_loop)
soft_dh_mutator   = Mutators.Dihedral.DihedralMutator(nb_dhs, () -> (randn() * soft_dh_mutator.step_size),  0.025, π/200)
soft_cs_mutator   = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dhs, () -> (randn() * soft_cs_mutator.step_size), 0.007, π/200)
soft_br_mutator   = Mutators.Blockrot.BlockrotMutator(metadata.blocks, () -> (randn() * soft_br_mutator.step_size), 0.15, π/100, 0.05, 15, loop_closer)

function dh_cs_br_soft_sampler!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, soft_dh_mutator)
    cm = Mutators.Crankshaft.run!(st, soft_cs_mutator)
    bm = Mutators.Blockrot.run!(st, soft_br_mutator)
    sampler_minimizer.run!(st, sampler_minimizer)
    return Dict("d" => dm, "c" => cm, "b" => bm)
end

initial_minimizer.run!(state, initial_minimizer)
refine_driver   = Drivers.MonteCarlo.Driver(dh_cs_br_soft_sampler!, amber_cr_evaluator!, 300.0, 100,
    print_status_refnm, print_structure_best, print_structure_rfnm)
refine_driver.run!(state, refine_driver)