using ProtoSyn
using ProtoSyn.Peptides

# "Pseudo-backbone" form
for aa in ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    pose = Peptides.build(Peptides.grammar, [aa])
    Peptides.cap!(pose)
    ProtoSyn.pop_atom!(pose, pose.graph[1, 1, 2])
    N = length(pose.graph[1, 1].items)
    pose.graph[1, 1, 3].symbol = "C"; pose.graph[1, 1, 3].name = "CX"
    pose.graph[1, 1, 2].symbol = "H"; pose.graph[1, 1, 2].name = "H"
    pose.graph[1, 1, N].symbol = "N"; pose.graph[1, 1, N].name = "NX"
    pose.graph[1, 1, N-1].symbol = "O"; pose.graph[1, 1, N-1].name = "O"
    ProtoSyn.write(pose, "$(lowercase(Peptides.one_2_three[aa[1]])).pdb")
end