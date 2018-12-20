using ProtoSyn

state, metadata     = Common.load_from_pdb("hydrogen_bonding_test.pdb")
hb_groups           = Forcefield.CoarseGrain.compile_hb_groups(metadata.atoms, 1.0)
eH                  = Forcefield.CoarseGrain.evaluate!(hb_groups, state)
println("Hydrogen Bonding Energy: $eH")