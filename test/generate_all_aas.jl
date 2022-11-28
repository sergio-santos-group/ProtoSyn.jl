using ProtoSyn
pose = ProtoSyn.Peptides.build(ProtoSyn.Peptides.grammar, seq"ECDARGNQPMKFIHWTSLYV")
ProtoSyn.write(pose, "all_aas.pdb")