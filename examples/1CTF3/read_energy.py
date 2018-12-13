from matplotlib import pyplot as plt
import sys

steps, energy = [], []
with open(sys.argv[1], "r") as fin:
    for line in fin:
        elem = line.split()
        steps.append(int(elem[1]))
        energy.append(float(elem[4]))

plt.plot(steps, energy)
plt.show()
