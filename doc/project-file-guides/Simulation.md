# Simulation.py

## Purpose

`Simulation.py` defines the `simulation` class, which owns the main run loop, recorded datasets, summary output, and final graph generation.

## Constructor responsibilities

When a `simulation` object is created, it:

- stores the simulation name and input parameters
- normalises the box dimensions into a numeric 3D size
- builds a `ParticleArray` wrapper around the full particle list
- computes the box volume and total wall area
- creates the graph-output directory for the run

## Main run loop

The `run_simulation()` method is the core of the programme.

### Step 1: convert time units

Time settings are converted into seconds so the runtime loop has consistent numeric values.

### Step 2: create data buffers

The simulation allocates lists for all sampled datasets, including:

- positions
- velocities
- translational energies
- temperatures
- pressures
- collision counts
- pairwise distances
- nearest-neighbour distances
- angular momentum
- rotational energy
- vibrational energy

### Step 3: record the initial state

The simulation records the state at `t = 0` before advancing any particles.

### Step 4: advance epochs

For each epoch, the simulation:

1. updates particle positions
2. resolves collisions
3. samples data every `record_interval`

### Step 5: finalise

After the loop, it:

- syncs the vectorised array state back into the original particle objects
- prints angular-momentum summary statistics
- prints energy and temperature summary statistics
- calls `plot_graphs(...)` to write the final figures

## Plot generation

`plot_graphs(...)` creates one graph object per output figure and points each one at the current run's save directory.

## Important practical detail

This file stores many recorded datasets across time. The pairwise-distance dataset is especially expensive because its size scales approximately like $O(N^2)$ per sample. Large particle counts can therefore become memory-limited.

## When to edit this file

Edit `Simulation.py` when changing the run loop, sampling schedule, summary statistics, or which graphs are produced at the end of a run.