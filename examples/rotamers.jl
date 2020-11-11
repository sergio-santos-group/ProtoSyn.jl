using ProtoSyn
using ProtoSyn.Peptides.Rotamers

rl = Rotamers.load_dunbrack(ProtoSyn.resource_dir * "/Peptides/dunbrack_rotamers.lib")
println(rl)