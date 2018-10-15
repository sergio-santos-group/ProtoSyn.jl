from matplotlib import pyplot as plt
from math import ceil

data = []
with open("mov.txt", "r") as f:
    for line in f:
        data.append(ceil(float(line)/2))

plt.plot(data)