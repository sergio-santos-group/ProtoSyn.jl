"""
    plot_energy.py

Read an energy.out file (containing all the energy components of the BEST structures from ProtoSyn sampling/refinement stages) and plot the energy values.
Dotted lines represent the reference structure energy.
"""

from matplotlib import pyplot as plt

keys = ["eTotal", "eBond", "eAngle", "eDihedral", "eCoulomb", "eCoulomb14", "eLJ", "eLJ14", "eContact", "eSol", "eDihedralFBR"]
scale_factor = [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k']
data, axis, threshold, targets = {}, {}, {}, {}
steps = []
for plot_index, key in enumerate(keys, start = 1):
    data[key] = []
    axis[key] = plt.subplot(len(keys), 1, plot_index)
    threshold[key] = None
    targets[key] = None

threshold_set = False
read_every = 1
fin = open("out/energy.out", "r")
for line_index, line in enumerate(fin):
    if line.startswith("(TRGT)"):
        try:
            elem = [float(x) for x in line.split()[1:]]
        except:
            continue
        for index, key in enumerate(keys):
            targets[key] = elem[index + 1]
    elif line_index % read_every == 0 and line.startswith("(BEST)"):
        try:
            elem = [float(x) for x in line.split()[1:]]
        except:
            continue
        steps.append(elem[0])
        if threshold_set == False:
            for index, key in enumerate(keys):
                threshold[key] = elem[index + 1] + abs(elem[index + 1]) * scale_factor[index]
            threshold_set = True
        for index, key in enumerate(keys):
            if elem[index + 1] >= threshold[key]:
                data[key].append(threshold[key])
            else:
                data[key].append(elem[index + 1])
fin.close()

for color_index, key in enumerate(keys):
    axis[key].axhline(targets[key], color = colors[color_index], linewidth = 2, linestyle = "--", alpha = 0.35)
    axis[key].scatter(steps, data[key], label = key, color = colors[color_index])
    axis[key].plot(steps, data[key], color = colors[color_index], label = None, linewidth = 0.75, alpha = 0.75)
    axis[key].legend(loc = 1)
plt.subplots_adjust(hspace = 0.1)
plt.show()