# Phys389 Particle Simulation

This repository contains a particle simulation written in Python for Phys389. The main workflow is interactive: you start the interface script, choose an existing initialisation file or enter parameters manually, the simulation runs, and plots are written to the plotted-graphs output folder.

## What the programme does

The simulation models particles moving inside a 3D box, handling:

- translational motion
- particle-particle collisions
- particle-wall collisions
- temperature, pressure, and energy diagnostics
- optional rotational and vibrational energy bookkeeping
- graph generation for the recorded time history

The code supports multiple particle species in one run, each with its own mass, radius, temperature, position distribution, and velocity distribution.

## Repository areas that matter for normal use

- `project/Interface.py`: interactive entry point
- `project/Initialiser.py`: reads parameters from file or user input
- `project/Distributions.py`: generates initial positions, velocities, and angular momenta
- `project/Particles.py`: particle data model and vectorised particle-array operations
- `project/Simulation.py`: main simulation loop and data recording
- `project/Graphs.py`: plotting and graph saving
- `project/initialisation_files/`: saved input files
- `project/plotted_graphs/`: generated graph output

Additional detailed notes for each main module are in `doc/project-file-guides/`.

## Requirements

The project declares these Python dependencies in `pyproject.toml`:

- `numpy`
- `scipy`
- `matplotlib`
- `sympy`

Python 3.8 or newer is required.

## Setup

From the repository root, create or activate a Python environment and install the dependencies you need. A typical pip-based setup is:

```powershell
python -m pip install numpy scipy matplotlib sympy
```

## Running the programme

The current import layout is most reliable when run from inside the `project/` directory.

```powershell
cd project
python Interface.py
```

When the programme starts, it asks whether to use file-based initialisation.

### Option 1: Use an existing initialisation file

1. Start the programme.
2. Enter `Y` when asked `Use file initialisation? (Y/N):`
3. Enter the name of a file from `project/initialisation_files/`, for example:

```text
water.txt
```

The simulation name will be derived from the input filename, and graphs will be saved into a matching subfolder under `project/plotted_graphs/`.

### Option 2: Enter parameters manually

If you choose not to load a file, the programme will prompt for:

1. box dimensions and length units
2. total simulation time, timestep, and graph-sampling interval
3. one or more particle definitions
4. an optional random seed
5. whether to save the custom initialisation to a file

If you skip saving the manual configuration, the programme will ask for a simulation name before running.

## Initialisation file format

Saved initialisation files are plain text files with comments and comma-separated values. In order, they contain:

1. box parameters: `Len_X,Units_X,Len_Y,Units_Y,Len_Z,Units_Z`
2. time parameters: `Time,Units_Time,dt,Units_dt,T_plt,Units_T_plt`
3. random seed
4. one or more particle lines

Particle lines use this format:

```text
type,N,m(Da),D,r,Temp,Pos_dist,Vel_dist,I_x,I_y,I_z,vib_modes,gradient_expr
```

Notes:

- `type`: particle label used in summaries and graphs
- `N`: number of particles of that type
- `m(Da)`: particle mass in daltons
- `D`: degrees of freedom
- `r`: collision radius
- `Temp`: initial temperature in kelvin
- `Pos_dist`: `1` lattice, `2` uniform, `3` custom gradient
- `Vel_dist`: `1` identical RMS, `2` uniform RMS range, `3` Maxwell-Boltzmann
- `I_x,I_y,I_z`: moments of inertia data used for rotational energy
- `vib_modes`: vibrational mode list
- `gradient_expr`: optional symbolic density expression for custom spatial sampling

## Outputs

During a run, the programme prints epoch progress and summary statistics to the terminal. At the end of the run it writes graphs to:

- `project/plotted_graphs/<simulation-name>/`

Typical outputs include:

- radial position distribution over time
- speed distribution over time
- translational energy distribution over time
- temperature distribution over time
- pressure history
- collision-count plots
- pairwise-distance plot
- nearest-neighbour-distance plot
- rotational, vibrational, and angular-momentum plots when those data are present

## Known practical limits

- Large particle counts can consume a lot of memory because pairwise-distance data are recorded across many sample times.
- Very long simulation names or very long graph titles can cause Windows path-length issues if file paths become too long.
- Because of the current imports in `project/Simulation.py`, launching from inside `project/` is the safer option.

## Where to read next

- `doc/project-file-guides/README.md` for the documentation index
- `doc/project-file-guides/Simulation.md` for the main runtime flow
- `doc/project-file-guides/Initialiser.md` for the input-file and manual-input format