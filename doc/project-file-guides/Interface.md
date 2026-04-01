# Interface.py

## Purpose

`Interface.py` is the interactive entry point for a simulation run. It ties together parameter loading, initial particle generation, simulation construction, and execution.

## Main flow

At startup, the file:

1. imports the other project modules
2. calls `INIT.Initialise()` to get box parameters, time parameters, particle definitions, a random seed, and optionally the source initialisation filename
3. seeds NumPy's random number generator
4. loops over each particle-type definition and creates the particle objects
5. chooses a simulation name
6. constructs `SIMU.simulation(...)`
7. runs the simulation and prints total runtime

## What it expects from the initialiser

The initialiser returns five values:

- `Box_Params`
- `Time_Params`
- `Sep_Particles`
- `seed`
- `init_file`

If `init_file` is present, `Interface.py` uses the filename as the simulation name. If not, it prompts the user to enter a simulation name manually.

## Particle construction

For each particle type, `Interface.py` asks `Distributions.py` to generate:

- starting positions
- starting velocities
- starting angular momenta

It then creates one `particle` object per particle and stores all of them in a single list for the simulation.

## Why this file matters

If you want to change how a run starts, how the simulation is named, or what gets built before the main loop begins, this is the first place to look.