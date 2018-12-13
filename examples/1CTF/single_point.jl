using ProtoSyn

state, metadata = Common.load_from_pdb("data/superposition_no_sc.pdb")
amber_topology  = Forcefield.Amber.load_from_json("data/1ctf_amber_top_no_sc.json")
contact_restraints  = Forcefield.Restraints.load_distance_restraints_from_file("data/contact_map_raptorx.txt", metadata, k = 10.0, threshold = 0.31)

function amber_cr_evaluator_nc!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    amber_e   = Forcefield.Amber.evaluate!(amber_topology.bonds, state, do_forces = do_forces)
    amber_e  += Forcefield.Amber.evaluate!(amber_topology.angles, state, do_forces = do_forces)
    amber_e  += Forcefield.Amber.evaluate!(amber_topology.atoms, state, do_forces = do_forces, cut_off = 1.2, eCoulomb_λ = 0.0)
    amber_e  += Forcefield.Amber.evaluate!(amber_topology.dihedralsCos, state, do_forces = do_forces)
    state.energy.comp["amber"] = amber_e
    contact_e = Forcefield.Restraints.evaluate!(contact_restraints, st, do_forces=do_forces)
    energy    = amber_e + contact_e
    st.energy.comp["other"] = contact_e
    st.energy.eTotal = energy
    return energy
end

amber_cr_evaluator_nc!(state, false)
Print.energy_by_component(state.energy)