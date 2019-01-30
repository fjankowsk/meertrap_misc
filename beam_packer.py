#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Compute beam packing for MeerTRAP.
#   2019 Fabian Jankowski
#

import numpy as np
import logging
import os.path

logger = logging.getLogger()

def load_data(filename):
    """
    Load the beam position data from file.
    """

    if not os.path.isfile(filename):
        raise RuntimeError("Input file does not exist: {filename}")
    
    dtype = [("x","float"), ("y","float")]

    data = np.genfromtxt(filename, delimiter="\t", dtype=dtype)

    return data


#
# MAIN
#

def main():
    infile = os.path.join("input", "134.0696_0.0_beam_pos.dat")
    data = load_data(infile)

    print(data)


if __name__ == "__main__":
    main()