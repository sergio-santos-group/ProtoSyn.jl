using ProtoSyn

state, metadata = Common.load_from_pdb("data/1ctf_native.pdb")
amber_topology  = Forcefield.Amber.load_from_json("data/1ctf_amber_top.json")
contact_restraints  = Forcefield.Restraints.load_distance_restraints_from_file("data/contact_map_raptorx.txt", metadata, k = 10.0, threshold = 0.31)

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