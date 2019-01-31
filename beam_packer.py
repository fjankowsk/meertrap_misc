#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Compute beam packing for MeerTRAP.
#   2019 Fabian Jankowski
#

import numpy as np
import logging
import os.path
import matplotlib.pyplot as plt
import sys
from timeit import default_timer as timer

logger = logging.getLogger()

def load_data(filename):
    """
    Load the beam position data from file.
    """

    if not os.path.isfile(filename):
        raise RuntimeError("Input file does not exist: {filename}")
    
    dtype = [("x","float"), ("y","float")]
    temp = np.genfromtxt(filename, delimiter="\t", dtype=dtype)

    dtype = [("nr","int"), ("x","float"), ("y","float"), ("group","int")]
    data = np.zeros(len(temp), dtype=dtype)

    data["nr"] = np.arange(len(data))
    for field in temp.dtype.names:
        data[field] = temp[field]

    return data


def plot_beam_centres(t_data):
    """
    Plot an overview of the beam centres.
    """

    data = np.copy(t_data)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(data["x"], data["y"],
               marker="x",
               color="black")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    fig.tight_layout()

    fig.savefig("beam_centres.pdf", bbox_inches="tight")
    fig.savefig("beam_centres.png", bbox_inches="tight",
                dpi=200)


def plot_beam_mapping(t_data):
    """
    Plot the beam mapping.
    """

    data = np.copy(t_data)

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    fig = plt.figure()
    ax = fig.add_subplot(111)

    i = 0

    for group in np.unique(data["group"]):
        color = colors[i % len(colors)]

        indiv = data[data["group"] == group]

        ax.scatter(indiv["x"], indiv["y"],
                marker="x",
                color=color,
                label=group)
        
        i += 1
    
    leg = ax.legend(loc="upper center", frameon=False, ncol=11, fontsize=9,
                    bbox_to_anchor=(0.5, -0.13),
                    columnspacing=0.2,
                    labelspacing=0.1)
    leg.set_zorder(10)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    fig.tight_layout()

    fig.savefig("beam_mapping.pdf", bbox_inches="tight")
    fig.savefig("beam_mapping.png", bbox_inches="tight",
                dpi=200)


def get_beam_mapping(t_data):
    """
    Map the beams to multicast addresses/compute nodes.
    """

    data = np.copy(t_data)

    nbeams = 396
    bunch = 6

    data = np.sort(data, order="x")

    if len(data) >= nbeams:
        data = data[0:nbeams]
        logger.info("Removed additional beams.")
    
    dtype = [("nr","int"), ("x","float"), ("y","float"), ("dist","float")]
    group = 0

    work = np.copy(data)

    while len(work) != 0:
        logger.debug("Length: {0}".format(len(work)))

        dist = np.sqrt((work["x"] - work["x"][0])**2 + (work["y"] - work["y"][0])**2)
        tot = np.zeros(len(work), dtype=dtype)

        for field in ["nr", "x", "y"]:
            tot[field] = work[field]

        tot["dist"] = dist

        tot = np.sort(tot, order="dist")

        # pick the closest six beams
        picked = tot[0:bunch]
        logger.debug("Group: {0}, beams: {1}".format(group, picked["nr"]))

        mask_work = np.isin(work["nr"], picked["nr"], invert=True)
        mask_data = np.isin(data["nr"], picked["nr"])
        work = work[mask_work]
        data["group"][mask_data] = group

        group += 1
    
    return data


#
# MAIN
#

def main():
    infile = os.path.join("input", "134.0696_0.0_beam_pos.dat")
    data = load_data(infile)

    plot_beam_centres(data)

    start = timer()
    mapped = get_beam_mapping(data)
    end = timer()

    print("Elapsed time: {0:.2f} ms".format(1000*(end - start)))

    for item in mapped:
        logger.info("Beam: {0}, group: {1}".format(item["nr"], item["group"]))

    plot_beam_mapping(mapped)

    plt.show()


if __name__ == "__main__":
    main()