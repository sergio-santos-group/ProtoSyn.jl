from protosyn.builder import grow, cap_chain

#Grow molecule from sequence
mol_sequence = "EFDVILKAAGANKVAVIKAVRGATGLGLKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK"
m = grow("mol", mol_sequence)
cap_chain(m)

def re_place(target, pos1, pos2):
    target.insert(pos2, target[pos1])
    target.pop(pos1)

for r in m.residues:
    # switch_pos(r.atoms, 0, 1)
    if r.index == 0: # Re-place H2 and H3 in correct order
        re_place(r.atoms, -1, 2)
        re_place(r.atoms, -1, 2)
m.renumber()

with open("output.pdb", "w") as file_out:
    file_out.write(m.as_pdb())