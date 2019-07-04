# ------------------ Energy Components -------------------------

λ_global = 1.0

# Component A: Amber
amber_top = Forcefield.Amber.load_from_json(input_amber_top)

# Component B: Coarse-grain solvation energy
solv_pairs = Forcefield.CoarseGrain.compile_solv_pairs(
    metadata.dihedrals,
    λ = 1.0 * λ_global)

# Component C: Coarse-grain hydrogen bonding network
hb_network = Forcefield.CoarseGrain.compile_hb_network(
    metadata.atoms,
    lib = Aux.read_JSON(hydrogen_bonding_library),
    λ = 0.4 * λ_global)

# Component D: Contacts distance restraints
contact_restraints = Forcefield.Restraints.load_distance_restraints_from_file(
    contacts_topology,
    metadata,
    λ = 800.0 * λ_global,
    threshold = 0.4,
    min_distance = 0.8)

# Component E: Dihedral restraints
dihedral_restraints = Forcefield.Restraints.lock_block_bb(metadata,
    fbw = 10.0,
    λ = 1e5)


# --------
# RFNM TEST

# Component B: Coarse-grain solvation energy (IN REFINEMENT)
solv_pairs_R = Forcefield.CoarseGrain.compile_solv_pairs(
    metadata.dihedrals,
    λ = 1.0)

# Component C: Coarse-grain hydrogen bonding network (IN REFINEMENT)
hb_network_R = Forcefield.CoarseGrain.compile_hb_network(
    metadata.atoms,
    lib = Aux.read_JSON(hydrogen_bonding_library),
    λ = 1.0)

# Component D: Contacts distance restraints (IN REFINEMENT)
contact_restraints_R = Forcefield.Restraints.load_distance_restraints_from_file(
    contacts_topology,
    metadata,
    λ = 800.0,
    threshold = 0.4,
    min_distance = 0.8)