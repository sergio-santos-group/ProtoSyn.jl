# -*- coding: utf-8 -*-
import argparse
import sys
import re

class Atom:
    """
    Holds information regarding one atom.
    """

    def __init__(self, id=-1, name='NaN', res_name='NaN', chain_name='A', res_id=-1, x=0.0, y=0.0, z=0.0, occ=1.0, tf=0.0, elem='X'):
        self.id         = int(id)
        self.name       = name
        self.res_name   = res_name
        self.chain_name = chain_name
        self.res_id     = int(res_id)
        self.x          = float(x)
        self.y          = float(y)
        self.z          = float(z)
        self.occ        = float(occ)
        self.tf         = float(tf)
        self.elem       = elem

    def as_pdb(self):
        """
        Print this atom information is PDB format.
        """
        return "%-6s%5d  %-4s% -3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n" % ("ATOM", self.id, self.name, self.res_name, self.chain_name,
            self.res_id, self.x, self.y, self.z, self.occ, self.tf, self.elem)


class Conect:
    """
    Holds information regarding the connectivity of an atom.
    """

    def __init__(self, id=-1, *conects):
        self.id         = int(id)
        self.conects    = map(int, filter(lambda c: c != '', conects))
    
    def as_pdb(self):
        """
        Print this connection information is PDB format.
        """
        return "%-6s%5d" % ("CONECT", self.id) + ''.join(["%5d" % m for m in self.conects]) + '\n'


class Protein:
    """
    Holds information about the protein object (atoms and connectivities).
    """

    def __init__(self, atoms=[], conects=[]):
        self.atoms      = atoms
        self.conects    = conects

    def as_pdb(self, title = None, out_file = None):
        """
        Print this protein information in PDB format
        to stdout, by default, an output file can be created by supplying `out_file`.
        Add a `title`.
        Usage:
            protein.as_pdb()
            protein.as_pdb("Sample 1", "sample_1.pdb")
        """
        if title != None:
            pdb_file = "TITLE   %-70s\n" % (title)
        else:
            pdb_file = ""
        for atom in self.atoms:
            pdb_file += atom.as_pdb()
        for conect in self.conects:
            pdb_file += conect.as_pdb()
        pdb_file += "END"
        if out_file == None:
            print pdb_file
        else:
            with open(out_file, "w") as file_out:
                file_out.write(pdb_file)
            print " | Wrote output PDB to %s file" % (out_file)
        
    
    def load_from_pdb(self, file_name):
        """
        Load atom and conects information from a `file_name` PDB file.
        Usage:
            protein.load_from_pdb("input_file.pdb")
        """
        pattern = re.compile('ATOM\s+([0-9]+)\s+([A-Z]*[0-9]*)\s+([A-Z]{3})\s+([A-Z]{1})\s*([0-9]+)\s+([-]?[0-9]*\.?[0-9]*)\s+([-]?[0-9]*\.?[0-9]*)\s+([-]?[0-9]*\.?[0-9]*)\s+([0-9]*\.?[0-9]*)\s+([0-9]*\.?[0-9]*)\s{6}\w?\s+(\w)', flags=re.S)
        with open(file_name, 'r') as file_in:
            results = pattern.findall(file_in.read())
        self.atoms = [Atom(*atom_info) for atom_info in results]
        pattern = re.compile('CONECT\s+(\d*)\s*(\d*)\s*(\d*)\s*(\d*)', flags=re.S)
        with open(file_name, 'r') as file_in:
            results = pattern.findall(file_in.read())
        self.conects = [Conect(*conect_info) for conect_info in results]
        print " | Loaded PDB object from %s file" % (file_name)


    def uniform_names_as_elem(self):
        """
        Sets all atom.name variables to the corresponding atom.elem (except CAs).
        Usage:
            protein.uniform_names_as_elem()
        """
        for atom in self.atoms:
            if atom.name != 'CA':
                atom.name = atom.elem
        print " | Set atom names as the corresponding elements"

    def regroup_hydrogens(self):
        """
        If there are is only 1 hydrogen per residue, reorder atoms so that (in the same order), there's 1 hydrogen in each residue.
        (Following Caterpillar model, ignoring Proline residues)
        Usage:
            protein.regroup_hydrogens()
        """
        hydrogens  = filter(lambda atom: atom.elem == 'H', self.atoms)
        n_residues = len(set([atom.res_id for atom in filter(lambda a: a.res_name != 'PRO', self.atoms)]))
        if len(hydrogens) != n_residues:
            print "ERROR: The loaded PDB file containts %d hydrogens and %d residues. Make sure theres only 1 hydrogen for each residue." % (len(hydrogens), n_residues)
            exit(1)
        self.atoms = filter(lambda atom: atom.elem != 'H', self.atoms)
        for index, atom in enumerate(self.atoms):
            if atom.elem == 'N' and atom.res_name != 'PRO':
                self.atoms.insert(index + 1, hydrogens[0])
                hydrogens.pop(0)
        print " | Regrouped hydrogen atoms 1 per residue"
        
    def renumber_atoms(self):
        """
        Renumber atoms (and corresponding connectivities) in the present order.
        Usage:
            protein.renumber_atoms()
        """
        conv = {}
        for index, atom in enumerate(self.atoms, start = 1):
            conv[atom.id] = index
            atom.id = index
        for atom in self.conects:
            atom.id = conv[atom.id]
            atom.conects = [conv[conect] for conect in atom.conects]
        print " | Renumbered atoms"

    def renumber_residues(self):
        """
        Renumber residues in the present order.
        Usage:
            protein.renumber_residues()
        """
        cur_res = -1
        count   = 0
        for atom in self.atoms:
            if atom.res_id != cur_res:
                cur_res = atom.res_id
                count += 1
            atom.res_id = count
        print " | Renumbered residues"


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Manipulate PDB files')
    parser.add_argument('-f', '--file', type=str, help='Input PDB file to be manipulated', required = True)
    parser.add_argument('-t', '--title', type=str, help='Add/Change PDB file title')
    parser.add_argument('-o', '--out', type=str, help='Output file: export output PDB to file')
    parser.add_argument('-rn', help='Rename: Change all atom names to the corresponding element (except CAs)', action = 'store_true')
    parser.add_argument('-gh', help='Group hydrogens (coarse-grain mode only): Distribute, in the present order, all hydrogens, 1 per residue;', action = 'store_true')
    parser.add_argument('-ra', help='Renumber atoms', action = 'store_true')
    parser.add_argument('-rr', help='Renumber residues', action = 'store_true')
    args = parser.parse_args()

    protein = Protein()
    protein.load_from_pdb(args.file)
    if args.rn: protein.uniform_names_as_elem()
    if args.gh: protein.regroup_hydrogens()
    if args.ra: protein.renumber_atoms()
    if args.rr: protein.renumber_residues()
    protein.as_pdb(args.title, args.out)