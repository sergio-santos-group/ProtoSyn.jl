using ProtoSyn
using Printf

state, metadata    = Common.load_from_pdb("data/1ctf_native_fixed.pdb")
amber_topology     = Forcefield.Amber.load_from_json("data/1ctf_amber_top.json")
contact_restraints = Forcefield.Restraints.load_distance_restraints_from_file("data/contact_map_raptorx.txt", metadata, k = 10.0, threshold = 0.31)
const out_pdb      = open("out/minimization.pdb", "w")
const out_ene      = open("out/minimization.dat", "w")

function amber_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    return Forcefield.Amber.evaluate!(amber_topology, st, cut_off = 1.2, do_forces=do_forces)
end

print_status = @Common.callback 1 function _print_status(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.Driver, max_force::Float64, gamma::Float64, args...)
    Print.status(@sprintf("(%5s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Max Force: %10.3e | Gamma: %10.3e\n", "I_SD", step, st.energy.eTotal,
        st.energy.comp["amber"], st.energy.comp["other"], max_force, gamma))
end

print_structure = @Common.callback 1 function _print_structure(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    Print.as_pdb(out_pdb, st, metadata, step = step)
    flush(out_pdb)
end

print_energy = @Common.callback 1 function _print_energy(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    println(out_ene, st.energy.eTotal)
    flush(out_ene)
end

initial_minimizer = Drivers.SteepestDescent.Driver(amber_evaluator!, 2000, 1000.0, 0.01, print_status, print_structure, print_energy)
initial_minimizer.run!(state, initial_minimizer)

function amber_cr_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    amber_e   = Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
    contact_e = Forcefield.Restraints.evaluate!(contact_restraints, st, do_forces=do_forces)
    energy    = amber_e + contact_e
    st.energy.eTotal = energy
    return energy
end

amber_cr_evaluator!(state, false)
Print.energy_by_component(state.energy)