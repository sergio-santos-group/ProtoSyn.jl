# -*- coding: utf-8 -*-
"""
     Created   José Pereira     April      2019
Last updated   José Pereira     April      2019
"""
from sys import stdout
import argparse
import re
"""
__________________________________________________________________________________________________________
clean_topol.py reads a "topol.top" file and performs the following functions:
- Turns N3 atoms to N.
- Turns O2 atoms to O and only accounts for 1 of them.
- Removes all "HP" and "H1".
- Only accounts for 1 H in each residue, removing all others.
- Renumbers all atoms to account for removed atoms, fixing the respective bonds, angles and dihedrals.
__________________________________________________________________________________________________________
"""

class Atom:

    def __init__(self, id = -1, name = 'NaN', res_id = -1, res_name = 'NaN', elem = 'NaN', cgnr = -1, charge = -1.0, mass = -1.0):
        self.id          = int(id)
        self.original_id = int(id)
        self.name        = name
        self.res_id      = int(res_id)
        self.res_name    = res_name
        self.elem        = elem
        self.cgnr        = int(cgnr)
        self.charge      = float(charge)
        self.mass        = float(mass)

    def __str__(self):
        pattern = "%6d %10s %6d %6s %6s %6d %10.4f %10.2f\n"
        return pattern % (self.id, self.name, self.res_id, self.res_name, self.elem, self.cgnr, self.charge, self.mass)


class Topology:

    def __init__(self, file_name):
        self.file_name = file_name
        self.atoms     = self.load_topology(file_name)

    def __getitem__(self, index):
        for atom in self.atoms:
            if atom.original_id == index:
                return atom.id

    def load_topology(self, file_name):
        pattern = re.compile('([0-9]+)\s+([A-Z]+[0-9]*)\s+([0-9]+)\s+([A-Z]+)\s+([A-Z]+[0-9]*)\s+([0-9]+)\s+([-]?[0-9]*\.?[0-9]*)\s+([-]?[0-9]*\.?[0-9]*)', flags=re.S)
        with open(file_name, 'r') as file_in:
            results = pattern.findall(file_in.read())
        return [Atom(*atom_info) for atom_info in results]

    def rename(self, target, new, attribute = 'name'):
        results = [(atom_index, atom) for atom_index, atom in enumerate(self.atoms) if getattr(atom, attribute) == target]
        for (atom_index, atom) in results:
            atom.name = new
            atom.elem = new
    
    def remove(self, target, attribute = 'name', leave = 0):
        results = [(atom_index, atom) for atom_index, atom in enumerate(self.atoms) if getattr(atom, attribute) == target]
        for (atom_index, atom) in reversed(results[:len(results)-leave]):
            self.atoms.pop(atom_index)

    def remove_per_residue(self, target, attribute = 'name', leave = 0):
        residue_list = set([atom.res_id for atom in self.atoms])
        for residue in residue_list:
            results = [(atom_index, atom) for atom_index, atom in enumerate(self.atoms) if getattr(atom, 'res_id') == residue]
            results = [(atom_index, atom) for atom_index, atom in results if getattr(atom, attribute) == target]
            for (atom_index, atom) in reversed(results[:len(results)-leave]):
                self.atoms.pop(atom_index)

    def renumber(self, start = 1):
        for (index, atom) in enumerate(self.atoms):
            atom.id = index + start

    def update_numbering(self, elem):
        s = ' '
        for atom_index in [int(e) for e in elem[:-1]]:
            if self[atom_index] == None: return None
            s += '%6s' % str(self[atom_index])
        return s + '%6s\n' % elem[-1]

    def print_from(self, file_name, output_stream = stdout):
        with open(file_name, 'r') as fin:
            for line in fin:
                elem = line.split()
                try:
                    int(elem[0])
                    if not re.compile('\s+[0-9]+\s+[A-Z]+[0-9]*\s+[0-9]+\s+[A-Z]+\s+[A-Z]+[0-9]*\s+[0-9]+\s+[-]?[0-9]*\.?[0-9]*\s+[-]?[0-9]*\.?[0-9]*').match(line):
                        s = self.update_numbering(elem)
                        if s is not None: output_stream.write(s)
                except:
                    output_stream.write(line)
                    if re.compile(';\sresidue\s+([0-9]+)').match(line):
                        residue = int(re.compile(';\sresidue\s+([0-9]+)').findall(line)[0])
                        for atom in [a for a in self.atoms if getattr(a, 'res_id') == residue]:
                            output_stream.write(atom.__str__())
                    continue
    
    def clean(self, output_stream = stdout):
        self.remove('H1')
        self.remove('HP')
        self.remove('O2', leave = 1)
        self.remove_per_residue('H', leave = 1)
        self.rename('O2', 'O')
        self.rename('N3', 'N')
        self.rename('H', 'H')
        self.renumber()
        self.print_from(self.file_name, output_stream)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Manipulate Topology files')
    parser.add_argument('-f', '--file', type=str, help='Input TOPOLOGY file to be manipulated', required = True)
    parser.add_argument('-o', '--out', type=str, help='Output TOPOLOGY destidation (Default: stdout)', default = stdout)
    args = parser.parse_args()

    topol = Topology(args.file)
    topol.clean(open(args.out, 'w'))
