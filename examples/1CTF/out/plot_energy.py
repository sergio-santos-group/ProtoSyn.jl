from matplotlib import pyplot as plt

keys = ["eTotal", "eBond", "eAngle", "eDihedral", "eCoulomb", "eCoulomb14", "eLJ", "eLJ14", "other"]
scale_factor = [-0.1, 1.0, 1.0, 0.5, 0.5, 0.5, 0.0, 0.0, -0.5]
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k']
data, axis, threshold = {}, {}, {}
for plot_index, key in enumerate(keys, start = 1):
    data[key] = []
    axis[key] = plt.subplot(len(keys), 1, plot_index)
    threshold[key] = None

threshold_set = False
fin = open("out/energy.out", "r")
for line in fin:
    try:
        elem = [float(x) for x in line.split()]
    except:
        continue
    if threshold_set == False:
        for index, key in enumerate(keys):
            threshold[key] = elem[index] + abs(elem[index]) * scale_factor[index]
        threshold_set = True
    for index, key in enumerate(keys):
        if elem[index] >= threshold[key]:
            data[key].append(threshold[key])
        else:
            data[key].append(elem[index])
fin.close()

for color_index, key in enumerate(keys):
    axis[key].plot(data[key], label = key, color = colors[color_index])
    axis[key].legend(loc = 1)
plt.subplots_adjust(hspace = 0.1)
plt.show()