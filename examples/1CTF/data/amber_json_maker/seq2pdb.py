from protosyn.molecule import Molecule
from protosyn.builder import grow, cap_chain
import StringIO

"""
    python seq2pdb.py > mol.pdb

seq2pdb.py script is responsible for producing straight mol.pdb files from the sequence.
"""

#Grow molecule from sequence
mol_sequence = "HRQALGERLYPRVQAMQPAFASKITGMLLELSPAQLLLLLASEDSLRARVDEAMELIIAHG"
m = grow("mol", mol_sequence)

#Print output as PDB
print m.as_pdb(include_bonds=True)
