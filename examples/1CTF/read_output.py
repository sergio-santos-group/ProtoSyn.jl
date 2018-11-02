from matplotlib import pyplot as plt

steps, energy, pmut, ar = [], [], [], []
t = None
with open("out/data.out") as f:
    for line in f:
        elem = line.split()
        steps.append(int(elem[2]))
        energy.append(float(elem[5]))
        pmut.append(float(elem[10]))
        ar.append(float(elem[18]))
        if t == None: t = float(elem[21])

grid = plt.GridSpec(3, 1, wspace=0.4, hspace=0.5)
ax1 = plt.subplot(grid[0, 0], title = "Energy (T = %4.2f)" % t)
ax2 = plt.subplot(grid[1, 0], title = "P_muts")
ax3 = plt.subplot(grid[2, 0], title = "Acceptance ratio")
ax1.plot(steps, energy)
ax2.plot(steps, pmut)
ax3.plot(steps, ar)
plt.show()