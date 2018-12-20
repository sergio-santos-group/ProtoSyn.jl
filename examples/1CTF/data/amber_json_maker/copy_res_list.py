"""
    copy_res_list.py template.pdb target.pdb

Copy a list of residues from template.pdb and paste it on target.pdb.
Both structures must have the same number of atoms.
"""

import sys

old_res = None
res_list = []
with open(sys.argv[1], "r") as fin:
    for line in fin:
        if line.startswith("ATOM"):
            elem = line.split()
            if elem[5] != old_res:
                old_res = elem[5]
                res_list.append(elem[3])

old_res = None
index = -1
with open(sys.argv[2], "r") as fin:
    for line in fin:
        if line.startswith("ATOM"):
            elem = line.split()
            if elem[5] != old_res:
                index += 1
                old_res = elem[5]
            print("ATOM {:6d}  {:<3s} {:3s} A {:3d} {:11.3f} {:7.3f} {:7.3f}  1.00 {:5.2f}           {:3s}".format(int(elem[1]), elem[2][:3], res_list[index], index,
                float(elem[6]), float(elem[7]), float(elem[8]), float(elem[10]), elem[11]))
        else:
            print(line[:-1])