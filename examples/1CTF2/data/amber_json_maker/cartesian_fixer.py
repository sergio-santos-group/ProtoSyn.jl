hs = []
with open("teste.xyz", "r") as fin:
    with open("teste_out.xyz", "w") as fout:
        for line in fin:
            if line.startswith("H"):
                hs.append(line)
            else:
                fout.write(line)
                if line.startswith("N"):
                    fout.write(hs[0])
                    hs.pop(0)