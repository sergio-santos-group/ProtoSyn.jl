from protosyn.dihedral import DihedralType as DType
from protosyn.molecule import Molecule
from collections import OrderedDict
from protosyn.builder import grow, cap_chain
import StringIO
import json

"""
    CONV.PY

    python conv.py > mol.json

    conv.py script is responsible for producing mol.json and mol.pdb files, two requesites for ILSRR-BCS
algorithm in Catulia. This script uses the Protosyn environment to create a molcule object contianing
the residues and dihedrals definitions. The necessary values are then printed to a JSON file that can
be parsed by Catulia. mol.pdb is a PDB model that is adapted by Catulia to print the produced structures.
"""

#Grow molecule from sequence
mol_sequence = "EFDVILKAAGANKVAVIKAVRGATGLGLKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK"
ss           = "CEEEEEEECCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCEEEEEEECHHHHHHHHHHHHHHCCEEEEC"
m = grow("mol", mol_sequence)
print len(m.get_coordinates())

m2 = Molecule.LoadFromFile("1ctf_native.pdb")
cap_chain(m)
print len(m2.get_coordinates())
m.set_coordinates(m2.get_coordinates())



def re_place(target, pos1, pos2):
    target.insert(pos2, target[pos1])
    target.pop(pos1)

for r in m.residues:
    # switch_pos(r.atoms, 0, 1)
    if r.index == 0: # Re-place H2 and H3 in correct order
        re_place(r.atoms, -1, 2)
        re_place(r.atoms, -1, 2)
m.renumber()

#Reload molecule to update dihedrals and connectivities
output = StringIO.StringIO()
output.write(m.as_pdb())
m = Molecule.LoadFromFileHandle(output)
output.close()

#Renumerate atoms
for r in m.residues:
    l = len(r.atoms)
    for dtype,dihd in r.dihedrals.iteritems():
        if dtype == DType.PHI:
            dihd.movable = range(3, l)
        elif dtype == DType.PSI:
            dihd.movable = range(l - 2, l)

#Create residue dictionary
residues = []
for r in m.residues:
    if r.next == None:
        _next = None
    else:
        _next = r.next.index + 1
    residues.append({
        'n': r.index + 1,
        'type': r.letter,
        'atoms': [x.index + 1 for x in r.atoms],
        'next': _next,
        'ss': ss[r.index]
    })

#Create dihedrals dictionary
dihedrals = []
for r in m.residues:
    keys = r.dihedrals.keys()
    for dtype,dihd in r.dihedrals.iteritems():
        extra = []
        if dtype == DType.OMEGA:
            continue

        movable = [index for index in dihd.movable]
        movable = [r.atoms[x] for x in movable]
        dihedrals.append({
            'a1': r.get_atom_by_name(dihd.a1).index + 1,
            'a2': r.get_atom_by_name(dihd.a2).index + 1,
            'a3': r.get_atom_by_name(dihd.a3).index + 1,
            'a4': r.get_atom_by_name(dihd.a4).index + 1,
            'parent' : r.index + 1,
            'movable' : [x.index + 1 for x in movable],
            'type' : dtype.name
        })

#Print output as JSON
with open('mol.pdb', 'w') as fout: fout.write(m.as_pdb(include_bonds=True))
print json.dumps(dict(residues=residues, dihedrals=dihedrals), indent=2)
