atom_count = 1
h_count    = 0
o2_count   = 0
conv = {}
with open("topol.top", "r") as fin:
    for line in fin:
        elem = line.split()
        try:
            int(elem[0])
        except:
            print line[:-1]
            h_count = 0
            continue
        if len(elem) == 11:
            if elem[1] in ["HP", "H1"]:
                    continue
            elif elem[1] == "N3":
                print "%6d %10s %6d %6s %6s %6d %10.4f %10.2f   ; qtot %6.4f" % (atom_count, "N", int(elem[2]), elem[3], "N", atom_count, float(elem[6]), float(elem[7]), float(elem[10]))
                conv[int(elem[0])] = atom_count
            elif elem[1] == "O2":
                o2_count += 1
                if o2_count > 1:
                    continue
                print "%6d %10s %6d %6s %6s %6d %10.4f %10.2f   ; qtot %6.4f" % (atom_count, "O", int(elem[2]), elem[3], "O", atom_count, float(elem[6]), float(elem[7]), float(elem[10]))
                conv[int(elem[0])] = atom_count
            else:
                if elem[1] == "H":
                    h_count += 1
                    if h_count > 1:
                        continue
                print "%6d %10s %6d %6s %6s %6d %10.4f %10.2f   ; qtot %6.4f" % (atom_count, elem[1], int(elem[2]), elem[3], elem[4], atom_count, float(elem[6]), float(elem[7]), float(elem[10]))
                conv[int(elem[0])] = atom_count
            atom_count += 1
        elif len(elem) == 3:
            try:
                print "%5d %5d %5d" % (conv[int(elem[0])], conv[int(elem[1])], int(elem[2]))
            except:
                continue
        elif len(elem) == 4:
            try:
                print "%5d %5d %5d %5d" % (conv[int(elem[0])], conv[int(elem[1])], conv[int(elem[2])], int(elem[3]))
            except:
                continue
        elif len(elem) == 5:
            try:
                print "%5d %5d %5d %5d %5d" % (conv[int(elem[0])], conv[int(elem[1])], conv[int(elem[2])], conv[int(elem[3])], int(elem[4]))
            except:
                continue