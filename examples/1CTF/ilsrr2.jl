using ProtoSyn

state, metadata = Common.load_from_pdb("data/1ctf_good_rmsd.pdb")

amber_topology = Forcefield.Amber.load_from_json("data/1ctf_amber_top.json")

nb_dihedrals       = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dihedrals   = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dihedrals) # For crankshaft movements
dihedral_mutator   = Mutators.Dihedral.DihedralMutator(nb_dihedrals, () -> (randn() * dihedral_mutator.step_size), 0.1, 0.1)
crankshaft_mutator = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dihedrals, () -> (randn() * crankshaft_mutator.step_size), 0.1, 0.1)
function my_sampler!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, dihedral_mutator)
    cm = Mutators.Crankshaft.run!(st, crankshaft_mutator)

    return Dict("d" => dm, "c" => cm)
end

inner_cycle_driver = Drivers.MonteCarlo.Driver(my_sampler!, my_evaluator!, temperature = 300.0, n_steps = 1)
ilsrr_driver = Drivers.ILSRR.Driver(inner_cycle_driver, my_evaluator!, my_sampler!)
ilsrr_driver.run!(state, ilsrr_driver)
#Same as: Drivers.ILSRR.run!(state, ilsrr_driver)