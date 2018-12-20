"""
    renumber_pdb.py input.pdb

Read an input.pdb and restart residue numbering from 1.
"""

from __future__ import print_function
import sys
index = 0
old_res = None
with open(sys.argv[1], "r") as fin:
    for line in fin:
        if line.startswith("ATOM"):
            elem = line.split()
            if elem[4] != old_res:
                index += 1
                old_res = elem[4]
            print("ATOM {:6d}  {:<3s} {:3s} A {:3d} {:11.3f} {:7.3f} {:7.3f}  1.00  0.00          {:3s}".format(int(elem[1]), elem[2][:3], elem[3], index,
                float(elem[5]), float(elem[6]), float(elem[7]), elem[10]))
        else:
            print(line[:-1])