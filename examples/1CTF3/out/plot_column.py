from matplotlib import pyplot as plt

data = []
with open("out/data.out", "r") as fin:
    for line in fin:
        elem = line.split()
        e = float(elem[5])
        if e <= 0:
            data.append(e)

plt.plot(data)
plt.show()