from matplotlib import pyplot as plt

steps, energy, max_force, gamma = [], [], [], []
with open("out/data.out") as f:
    for line in f:
        elem = line.split()
        if elem[1].startswith("SD"):
            steps.append(int(elem[2]))
            energy.append(float(elem[5]))
            max_force.append(float(elem[9]))
            gamma.append(float(elem[12]))

grid = plt.GridSpec(3, 1, wspace=0.4, hspace=0.5)
ax1 = plt.subplot(grid[0, 0], title = "Energy")
ax2 = plt.subplot(grid[1, 0], title = "Max force")
ax3 = plt.subplot(grid[2, 0], title = "Gamma")
ax1.plot(steps, energy)
ax2.plot(steps, max_force)
ax3.plot(steps, gamma)
ax3.axhline(0.2, c = "#3f7d20", linewidth = 3, alpha = 0.3, zorder = 1)
plt.show()