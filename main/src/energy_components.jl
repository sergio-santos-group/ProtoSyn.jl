# ------------------ Energy Components -------------------------

# Component A: Amber
amber_top = Forcefield.Amber.load_from_json(input_amber_top)

# Component B: Coarse-grain solvation energy
solv_pairs = Forcefield.CoarseGrain.compile_solv_pairs(
    metadata.dihedrals,
    λ = 1.0)

# Component C: Coarse-grain hydrogen bonding network
hb_network = Forcefield.CoarseGrain.compile_hb_network(
    metadata.atoms,
    lib = Aux.read_JSON(hydrogen_bonding_library),
    λ = 1_000_000.0)

# Component D: Contacts distance restraints
contact_restraints = Forcefield.Restraints.load_distance_restraints_from_file(
    contacts_topology,
    metadata,
    λ = 1000.0,
    threshold = 0.4,
    min_distance = 0.8)

# Component E: Dihedral restraints
dihedral_restraints = Forcefield.Restraints.lock_block_bb(metadata,
    fbw = 10.0,
    λ = 1e4)
