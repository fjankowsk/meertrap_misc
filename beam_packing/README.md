# Beam to compute node assignment optimisation code #

This directory contains code to determine the optimum assignment of coherent beams to compute nodes and IP multicast groups on the MeerTRAP cluster. The algorithm implemented optimises the spatial locality of beams on nodes, i.e. that neighbouring beams on the sky get processed on the same compute node. This allows for local (intra-node) multi-beam filtering, clustering, and sifting of single-pulse candidates.

There are implementations in `Mathematica` (by Sotiris Sanidas), `python` (FJ) and an unfinished version in `Golang` (FJ).

## Requirements ##

* Numpy
* Matplotlib
* Mathematica (for the initial implementation)
* [Golang](https://golang.org) (for the unfinished version)
