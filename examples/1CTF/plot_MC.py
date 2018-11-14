from matplotlib import pyplot as plt

steps, energy, amber, other, temp, ar = [], [], [], [], [], []
t = None
prev_steps = 0
init_energy, init_amber, init_other = None, None, None
energy_threshold = 10000
with open("out/data.out") as f:
    for line in f:
        elem = line.split()
        if len(elem) <= 1:
            continue
        elif elem[1].startswith("MC"):
            steps.append(prev_steps + int(elem[2]))
            if init_energy == None:
                init_energy = float(elem[5])
                init_amber  = float(elem[8])
                init_other  = float(elem[11])
            
            if float(elem[5]) <= init_energy + energy_threshold:
                energy.append(float(elem[5]))
            else:
                energy.append(init_energy + energy_threshold)

            if float(elem[8]) <= init_amber + energy_threshold:
                amber.append(float(elem[8]))
            else:
                amber.append(init_amber + energy_threshold)

            if float(elem[11]) <= init_other + energy_threshold:
                other.append(float(elem[11]))
            else:
                other.append(init_other + energy_threshold)

            temp.append(float(elem[27]))
            ar.append(float(elem[24]))
            if t == None:
                t = float(elem[27])
        elif elem[1].startswith("Perturbator"):
            prev_steps = steps[len(steps) - 1]

grid = plt.GridSpec(3, 1, wspace=0.4, hspace=0.5)
ax1 = plt.subplot(grid[0, 0], title = "Energy (T = %4.2f)" % t)
ax2 = plt.subplot(grid[1, 0], title = "Temperature")
ax3 = plt.subplot(grid[2, 0], title = "Acceptance ratio")
ax1.plot(steps, energy, label = "Total")
ax1.plot(steps, amber, label = "Amber")
ax1.plot(steps, other, label = "Other")
ax2.plot(steps, temp)
ax3.plot(steps, ar)
ax1.legend()
plt.show()