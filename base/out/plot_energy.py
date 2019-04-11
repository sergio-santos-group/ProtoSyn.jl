# -*- coding: utf-8 -*-
"""
     Created   José Pereira     February   2019
Last updated   José Pereira     March      2019
"""

from matplotlib import pyplot as plt
from operator import add
import numpy as np
import argparse

colors = ['blue', 'darkslategrey', 'sienna', 'darkorange', 'goldenrod', 'olivedrab', 'green', 'c', 'cadetblue', 'navy', 'mediumvioletred', 'crimson', 'pink']

def read_data_from_file(filename):
    keys = []
    data = {}
    with open(filename, "r") as fin:
        for line in fin:
            elem = line.split()
            if len(elem) > 0 and elem[0] == "Step":
                for key in elem[1:]:
                    keys.append(key)
                for key in keys:
                    data[key] = []
            elif len(elem) > 0 and elem[0] == "(BEST)":
                for index, value in enumerate(elem[2:]):
                    data[keys[index]].append(float(value))
    return keys, data


def plot_all(keys, data):
    global colors
    axis = {}
    for plot_index, key in enumerate(keys, start = 1):
        axis[key] = plt.subplot(len(keys), 1, plot_index)
        axis[key].scatter(range(0, len(data[key])), data[key], label = key, color = colors[plot_index])
        axis[key].plot(range(0, len(data[key])), data[key], label = None, linewidth = 0.75, alpha = 0.5, color = colors[plot_index])
        axis[key].legend(loc = 1, shadow = 1, bbox_to_anchor=(1.12, 1.0))

    plt.subplots_adjust(hspace = 0.1)
    plt.show()


def plot_grouped(keys, data):
    global colors
    amber_keys = ["eBond", "eAngle", "eDihedral", "eCoulomb", "eCoulomb14", "eLJ", "eLJ14"]
    cg_keys = ["eContacts", "eSol", "eH"]
    grouped_data = {"Amber": np.zeros(len(data[keys[0]])), "coarseGrain": np.zeros(len(data[keys[0]]))}
    for key in amber_keys:
        grouped_data["Amber"] += data[key]
    for key in cg_keys:
        grouped_data["coarseGrain"] += data[key]
    axis = {}
    for plot_index, key in enumerate(["Amber", "coarseGrain"], start = 1):
        axis[key] = plt.subplot(2, 1, plot_index)
        axis[key].scatter(range(0, len(grouped_data[key])), grouped_data[key], label = key, color = colors[plot_index])
        axis[key].plot(range(0, len(grouped_data[key])), grouped_data[key], label = None, linewidth = 0.75, alpha = 0.5, color = colors[plot_index])
        axis[key].legend(loc = 1, shadow = 1, bbox_to_anchor=(1.12, 1.0))
    
    plt.subplots_adjust(hspace = 0.1)
    plt.show()
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Plot energy contributions from ProtoSyn run')
    parser.add_argument('-f', '--file', type=str, help='Input file to be plotted', required = True)
    parser.add_argument('-g', '--group', help='Show only grouped contributions (Amber/CoarseGrain)', action='store_true')
    args = parser.parse_args()

    keys, data = read_data_from_file(args.file)

    if args.group:
        plot_grouped(keys, data)
    else:
        plot_all(keys, data)
