using ProtoSyn
using ProtoSyn.Peptides

# "Pseudo-backbone" form
for aa in ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "Q", "R", "S", "T", "V", "W", "Y"]
    pose = Peptides.build(Peptides.grammar, [aa])
    Peptides.cap!(pose)
    ProtoSyn.pop_atom!(pose, pose.graph[1, 1, 2])
    N = length(pose.graph[1, 1].items)
    pose.graph[1, 1, 3].symbol = "C"; pose.graph[1, 1, 3].name = "CX"
    pose.graph[1, 1, 2].symbol = "H"; pose.graph[1, 1, 2].name = "H"
    pose.graph[1, 1, N].symbol = "N"; pose.graph[1, 1, N].name = "NX"
    pose.graph[1, 1, N-1].symbol = "O"; pose.graph[1, 1, N-1].name = "O"
    ProtoSyn.write(pose, "$(lowercase(ProtoSyn.one_2_three[aa[1]])).pdb")
end

# Proline is a special case
for aa in ["P"]
    pose = Peptides.build(Peptides.grammar, [aa])
    CD = copy(pose.graph[1, 1, "CD"])
    CD_state = copy(pose.state[CD])
    Peptides.cap!(pose)
    ProtoSyn.pop_atom!(pose, pose.graph[1, 1, 2])
    N = length(pose.graph[1, 1].items)
    pose.graph[1, 1, 3].symbol = "C"; ProtoSyn.rename!(pose.graph[1, 1, 3], "CX")
    pose.graph[1, 1, 2].symbol = "H"; ProtoSyn.rename!(pose.graph[1, 1, 2], "H")
    pose.graph[1, 1, N].symbol = "N"; ProtoSyn.rename!(pose.graph[1, 1, N], "NX")
    pose.graph[1, 1, N-1].symbol = "O"; ProtoSyn.rename!(pose.graph[1, 1, N-1], "O")
    ProtoSyn.insert_atom_as_children!(pose, pose.graph[1, 1, "CG"], CD, CD_state)
    ProtoSyn.bond(CD, pose.graph[1, 1, "N"])
    ProtoSyn.pop_atom!(pose, pose.graph[1, 1].items[2])
    ProtoSyn.write(pose, "$(lowercase(ProtoSyn.one_2_three[aa[1]])).pdb")
end