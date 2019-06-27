using ProtoSyn

output_file = open("6btc_no_sc.pdb", "w")
state, metadata = Common.load_from_pdb("6btc_no_sc_native.pdb")
backbone = filter(x -> x.dtype <= Common.DIHEDRAL.omega, metadata.dihedrals)
Common.stretch_conformation!(state, backbone)
Print.as_pdb(output_file, state, metadata)
close(output_file)
