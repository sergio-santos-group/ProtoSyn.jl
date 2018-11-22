import sys

with open(sys.argv[1]) as fin:
    with open("stdr_contact_map.txt", "w") as fout:
        for line in fin:
            elem = line.split()
            if float(elem[4]) > 0:
                fout.write("%3s %3s %5s %5s %13s\n" % (elem[0], elem[2], elem[0] + "_" + elem[1], elem[2] + "_" + elem[3], elem[4]))
