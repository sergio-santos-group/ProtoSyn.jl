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
        data = {}
        with open(filename, "r") as file_in:
            l1 = file_in.readline()
            keys = filter(lambda x: len(x) > 1, l1.split())
            for key in keys:
                data[key] = []
            for line in file_in:
                ln = filter(lambda x: len(x) > 1, line[17:].split())
                for index, key in enumerate(keys):
                    point = float(ln[index])
                    if max_y != None and point > max_y:
                        data[key].append(max_y)
                    else:
                        data[key].append(point)
        return data


    def animate(i):
        data = parse_file("energy.log", 10000.0)
        ax1.clear()
        for key in data.keys():
            ax1.plot(data[key], label=key)
        ax1.legend()


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