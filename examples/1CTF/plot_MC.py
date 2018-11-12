from matplotlib import pyplot as plt

steps, energy, amber, other, pmut, ar = [], [], [], [], [], []
t = None
with open("out/data.out") as f:
    for line in f:
        elem = line.split()
        if elem[1].startswith("MC"):
            steps.append(int(elem[2]))
            energy.append(float(elem[5]))
            amber.append(float(elem[8]))
            other.append(float(elem[11]))
            pmut.append(float(elem[16]))
            ar.append(float(elem[24]))
            if t == None:
                t = float(elem[27])

grid = plt.GridSpec(3, 1, wspace=0.4, hspace=0.5)
ax1 = plt.subplot(grid[0, 0], title = "Energy (T = %4.2f)" % t)
ax2 = plt.subplot(grid[1, 0], title = "P_muts")
ax3 = plt.subplot(grid[2, 0], title = "Acceptance ratio")
ax1.plot(steps, energy, label = "Total")
ax1.plot(steps, amber, label = "Amber")
ax1.plot(steps, other, label = "Other")
ax2.plot(steps, pmut)
ax3.plot(steps, ar)
ax1.legend()
plt.show()