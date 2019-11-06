using ProtoSyn


# Creating a residue
residue = ProtoSyn.Residue(
    "BKB",
    [
        Atom(id=1,  name= "H", symbol="H",  x=0.3909, y=0.0724, z=0.0000),
        Atom(id=2,  name= "N", symbol="N",  x=0.3326, y=0.1548, z=0.0000),
        Atom(id=3,  name="CA", symbol="C",  x=0.3970, y=0.2846, z=0.0000),
        Atom(id=4,  name= "C", symbol="C",  x=0.5486, y=0.2705, z=0.0000),
        Atom(id=5,  name= "O", symbol="O",  x=0.6009, y=0.1593, z=0.0000),
        
    ],
    ConnectGraph(
        1 => [2],
        2 => [1, 3],
        3 => [2, 4],
        4 => [3, 5],
        5 => [4]
    )
)



# Creating a residue library (dictionary of resname=>residue)
lib = ProtoSyn.ResidueLib()
lib[residue.name] = residue
# or 
# ProtoSyn.ResidueLib(residue.name => residue)



#   A residue is a template for constructing molecules. An instance of
# residue within a molecule happens through instantiation of a LinkedResidue,
# having a residue as its source. Hence, a residue can be used multiple times
# (remember: it is a simple template).
#   A LinkedResidue has an id, a source, and can have multiple links (Link) to
# other (multiple) linked residues.
#   Example: two LinkedResidue(s) of type BKB (not linked to each other)
lr1 = ProtoSyn.LinkedResidue(id=1, source=residue)
lr2 = ProtoSyn.LinkedResidue(id=2, source=residue)


#   A molecule is a list of linked residues (LinkedResidue[]) and a
# list of links (Link[]) linking those residues. Nevertheless, links 
# between linked residues are not mandatory (may be changed in the future).
#   Example:
mol = ProtoSyn.Molecule(residues=[lr1, lr2])
