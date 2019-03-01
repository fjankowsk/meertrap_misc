#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Compute beam packing for MeerTRAP.
#   2019 Fabian Jankowski
#   Packing algorithm devised by Sotiris Sanidas 2019.
#

import numpy as np
import logging
import os.path
import matplotlib.pyplot as plt
import sys
from timeit import default_timer as timer

def load_data(filename):
    """
    Load the beam position data from file.

    Parameters
    ----------
    filename : str
        Name of the file from which to load the beam positions.
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


def plot_beam_packing(t_data):
    """
    Plot the beam packing.
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


def plot_packing_metric(t_data):
    """
    Visualize the packing metrics.
    """

    data = np.copy(t_data)

    data = np.sort(data, order="totdist")

    info_str = "min: {0:.2f}".format(np.min(data["totdist"]))
    info_str += ", median: {0:.2f}".format(np.median(data["totdist"]))
    info_str += ", max: {0:.2f}".format(np.max(data["totdist"]))
    info_str += ", tot: {0:.2f}".format(np.sum(data["totdist"]))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.hist(data["totdist"], bins=19,
            histtype="step",
            color="black",
            lw=2)

    ax.set_xlabel("Total distance (degrees)")
    ax.set_ylabel("Number")
    ax.grid(True)
    ax.set_title(info_str)

    fig.tight_layout()

    fig.savefig("packing_metric.pdf", bbox_inches="tight")
    fig.savefig("packing_metric.png", bbox_inches="tight",
                dpi=200)

    # cummulative
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.hist(data["totdist"], bins=19,
            histtype="step",
            color="black",
            lw=2,
            density=True, cumulative=True)

    ax.set_xlabel("Total distance (degrees)")
    ax.set_ylabel("CDF")
    ax.grid(True)
    ax.set_title(info_str)

    fig.tight_layout()

    fig.savefig("packing_metric_cum.pdf", bbox_inches="tight")
    fig.savefig("packing_metric_cum.png", bbox_inches="tight",
                dpi=200)


def get_beam_packing(t_data):
    """
    Map the beams to multicast addresses/compute nodes.
    """

    data = np.copy(t_data)

    nbeams = 396
    bunch = 6

    data = np.sort(data, order="x")

    logger = logging.getLogger()

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


def check_beam_packing(t_data):
    """
    Evaluate the beam packing using various metrics.
    """

    data = np.copy(t_data)

    data = np.sort(data, order="group")

    logger = logging.getLogger()

    dtype = [("group","int"), ("totdist","float")]
    info = np.zeros(np.max(data["group"]) + 1, dtype=dtype)

    for group in np.unique(data["group"]):
        indiv = data[data["group"] == group]

        indiv = np.sort(indiv, order="x")

        totdist = 0

        for nr in range(len(indiv) - 1):
            dist = np.sqrt((indiv["x"][nr+1:] - indiv["x"][nr])**2 + (indiv["y"][nr+1:] - indiv["y"][nr])**2)
            value = np.sum(dist)
            logger.debug("Group: {0}, nr: {1}, dist: {2}, val: {3}".format(group, nr, dist, value))

            totdist += value
            logger.debug("Total distance: {0}".format(totdist))
        
        info["group"][group] = group
        info["totdist"][group] = totdist

    return info


#
# MAIN
#

def main():
    logger = logging.getLogger()

    infile = os.path.join("input", "134.0696_90.0_beam_pos.dat")
    data = load_data(infile)

    plot_beam_centres(data)

    start = timer()
    packed = get_beam_packing(data)
    end = timer()

    print("Elapsed time: {0:.2f} ms".format(1000*(end - start)))

    for item in packed:
        logger.info("Beam: {0}, group: {1}".format(item["nr"], item["group"]))

    plot_beam_packing(packed)

    start = timer()
    metric = check_beam_packing(packed)
    end = timer()

    print("Elapsed time: {0:.2f} ms".format(1000*(end - start)))

    plot_packing_metric(metric)

    plt.show()


if __name__ == "__main__":
    main()