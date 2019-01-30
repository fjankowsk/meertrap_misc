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

logger = logging.getLogger()

def load_data(filename):
    """
    Load the beam position data from file.
    """

    if not os.path.isfile(filename):
        raise RuntimeError("Input file does not exist: {filename}")
    
    dtype = [("x","float"), ("y","float")]
    temp = np.genfromtxt(filename, delimiter="\t", dtype=dtype)

    dtype = [("nr","int"), ("x","float"), ("y","float"), ("picked","bool")]
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
        print("Length: {0}".format(len(work)))

        dist = np.sqrt((work["x"] - work["x"][0])**2 + (work["y"] - work["y"][0])**2)
        tot = np.zeros(len(work), dtype=dtype)

        for field in ["nr", "x", "y"]:
            tot[field] = work[field]

        tot["dist"] = dist

        tot = np.sort(tot, order="dist")

        # pick the closest six beams
        picked = tot[0:bunch]
        print("Group: {0}, beams: {1}".format(group, picked["nr"]))

        mask = np.isin(work["nr"], picked["nr"], invert=True)
        work = work[mask]
        group += 1


#
# MAIN
#

def main():
    infile = os.path.join("input", "134.0696_0.0_beam_pos.dat")
    data = load_data(infile)

    plot_beam_centres(data)
    get_beam_mapping(data)

    plt.show()


if __name__ == "__main__":
    main()