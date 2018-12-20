using ProtoSyn
using LinearAlgebra

state, metadata = Common.load_from_pdb("data/1ctf.pdb")
amber_topology  = Forcefield.Amber.load_from_json("data/1ctf_amber_top.json")
contact_restraints  = Forcefield.Restraints.load_distance_restraints_from_file("data/contact_map_raptorx.txt", metadata, k = 10.0, threshold = 0.31)

include("src/esol.jl")
const xyz_destination          = open("out/trajectory.pdb", "w")

function amber_cr_evaluator_nc!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    amber_e   = Forcefield.Amber.evaluate!(amber_topology.bonds, state, do_forces = do_forces)
    amber_e  += Forcefield.Amber.evaluate!(amber_topology.angles, state, do_forces = do_forces)
    amber_e  += Forcefield.Amber.evaluate!(amber_topology.atoms, state, do_forces = do_forces, cut_off = 1.2, eCoulomb_λ = 0.0)
    amber_e  += Forcefield.Amber.evaluate!(amber_topology.dihedralsCos, state, do_forces = do_forces)
    contact_e = Forcefield.Restraints.evaluate!(contact_restraints, st, do_forces=do_forces)
    sol_e     = evaluate!(solv_pairs, st)
    state.energy.comp["amber"] = amber_e
    st.energy.comp["other"]    = contact_e + sol_e
    st.energy.eTotal           = amber_e + contact_e + sol_e
    return amber_e + contact_e + sol_e
end

function amber_cr_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    amber_e   = Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
    contact_e = Forcefield.Restraints.evaluate!(contact_restraints, st, do_forces=do_forces)
    sol_e     = evaluate!(solv_pairs, st)
    state.energy.comp["amber"] = amber_e
    st.energy.comp["other"]    = contact_e + sol_e
    st.energy.eTotal           = amber_e + contact_e + sol_e
    return amber_e + contact_e + sol_e
end

amber_cr_evaluator!(state, false)
Print.energy_by_component(state.energy)
old_energy = deepcopy(state.energy)
Print.as_pdb(xyz_destination, state, metadata)

state, metadata = Common.load_from_pdb("data/1ctf_native.pdb")
amber_cr_evaluator!(state, false)
Print.energy_by_component(state.energy, old_energy)
Print.as_pdb(xyz_destination, state, metadata)
