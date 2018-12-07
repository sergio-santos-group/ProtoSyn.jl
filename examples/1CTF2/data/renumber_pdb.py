from __future__ import print_function
import sys
index = 1
conv = {}
with open(sys.argv[1], "r") as fin:
    for line in fin:
        if line.startswith("ATOM"):
            elem = line.split()
            print("ATOM {:6d}  {:<3s} {:3s} A {:3d} {:11.3f} {:7.3f} {:7.3f}  1.00  0.00 {:>11s}".format(index, elem[2], elem[3], int(elem[4]),
                float(elem[5]), float(elem[6]), float(elem[7]), elem[10]))
            conv[elem[1]] = index
            index += 1
        # elif line.startswith("CONECT"):
        #     elem = line.split()
        #     elem = elem[1:]
        #     print("\nCONECT", end="")
        #     for atom in [conv[atom] for atom in elem]:
        #         print("%4d", atom) 
        # else:
        #     print(line[:-1])