# Graphs.py

## Purpose

`Graphs.py` contains the plotting classes used at the end of a simulation run. Each graph class knows how to turn one recorded dataset into a Matplotlib figure and save it as a PNG.

## Base `graph` class

The base class stores shared metadata such as the title, axis labels, and sampling timestep. It also provides helper methods for:

- generating safe filenames for Windows
- choosing time-based x-axis values
- saving a figure into the current run's output directory

The save helper strips the simulation-name prefix from the graph title before writing the filename, which keeps saved paths shorter and avoids Windows path-length problems.

## Plot types in this file

The module includes graph classes for:

- position distributions
- velocity distributions
- translational energy distributions
- temperature distributions
- pressure histories
- collision counts per sample interval
- cumulative collision counts
- pairwise separation distances
- nearest-neighbour distances
- total angular momentum squared
- total rotational energy
- total vibrational energy
- combined energy-component histories

Most of the distribution-style plots are implemented as time-evolving histograms shown as heatmaps with a mean-value overlay.

## Data expectations

Each plot class expects the dataset shape produced by `Simulation.py`. For example:

- position and velocity plots expect one array per sampled epoch
- pressure plots expect ideal and impulse-based values for each particle
- pairwise-distance plots expect one 1D array of pair distances per sampled epoch

## Why this file matters

If graphs are missing, mislabeled, slow, or saved to the wrong place, this is the file to inspect first.

## When to edit this file

Edit `Graphs.py` when changing plot appearance, output filenames, graph types, axis conventions, or figure-saving behaviour.