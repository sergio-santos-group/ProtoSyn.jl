# -*- coding: utf-8 -*-
"""
     Created   José Pereira     June       2019
Last updated   José Pereira     June       2019
"""
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import argparse

style.use('seaborn-bright')

def main(static = False):

    def parse_file(filename, max_y = None):
        data = {"home_state": [], "best_state": []}
        current_home_state = 0.0
        current_best_state = 0.0
        rfnm_list = []
        rfnm_index = -2
        with open(filename, "r") as file_in:
            l1 = file_in.readline()
            keys = filter(lambda x: len(x) > 1, l1.split())
            for key in keys:
                data[key] = []
            for line in file_in:
                if line[:5] == "ILSRR":
                    ln = filter(lambda x: len(x) > 1, line[6:].split())
                    current_home_state = min(float(ln[0]), max_y)
                    current_best_state = min(float(ln[1]), max_y)
                else:
                    if line[:16].split()[0] == "RFNM":
                        rfnm_list.append(rfnm_index)
                    ln = filter(lambda x: len(x) > 1, line[17:].split())
                    data["home_state"].append(current_home_state)
                    data["best_state"].append(current_best_state)
                    for index, key in enumerate(keys):
                        data[key].append(min(float(ln[index]), max_y))
                    rfnm_index += 1
        return rfnm_list, data


    def animate(i):
        colors = ["cyan", "teal", "navy", "crimson", "magenta", "orange", "navy", "green"]
        dashes = ["-", "--", "--", "-", "-", "-", "-", "-"]
        alphas = [0.75, 1.0, 1.0, 0.75, 0.75, 0.75, 1.0, 1.0]
        rfnm_list, data = parse_file("energy.log", 6000.0)
        ax1.clear()
        for index, key in enumerate(data.keys()):
            ax1.plot(data[key][1:],
                dashes[index],
                label=key,
                color = colors[index],
                alpha = alphas[index])
        for index, key in enumerate(data.keys()):
            ax1.plot(len(data[key])-2, data[key][0],
                "o",
                color = colors[index])
        for rfnm_index in rfnm_list:
            ax1.axvspan(rfnm_index, rfnm_index + 1, facecolor = 'khaki', alpha = 0.6)
        ax1.legend(loc = 2, prop={'size': 6})


    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)

    animate(0)
    if static != True:
        ani = animation.FuncAnimation(fig, animate, interval=5000)
    plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Plot energy contributions from ProtoSyn run')
    parser.add_argument('-s', '--static',
        help="Don't animate the plot",
        action = 'store_true')
    args = parser.parse_args()

    print "\n   Written by José M. Pereira @ Universidade de Aveiro, 2019\n   Ploting data from 'energy.log' ..."
    main(static = args.static)