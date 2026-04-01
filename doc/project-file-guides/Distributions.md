# Distributions.py

## Purpose

`Distributions.py` generates the initial state of the simulation. It contains the code that samples particle positions, velocities, and angular momenta before the simulation begins.

## Position distributions

The file defines a common 3D coordinate-distribution base class and several concrete implementations.

### `Lattice_s`

Places particles on a regular lattice inside the box. It tries to factor the particle count into a grid that fits the box dimensions reasonably well.

### `Uniform_prob_density_s`

Samples positions uniformly throughout the box.

### `Gradient_prob_density_s`

Samples positions from a custom symbolic probability density function of `x`, `y`, and `z`.

This class:

- parses the symbolic expression using SymPy
- converts it to a NumPy-callable function
- evaluates the density on a coarse grid
- chooses cells with probability proportional to the density
- samples random positions inside the chosen cells

This is how the atmosphere-style gradient setups work.

## Velocity distributions

The file also defines a base class for 3D vector distributions and three implementations:

- identical RMS speed in random directions
- uniformly distributed RMS speed up to a temperature-based maximum
- Maxwell-Boltzmann sampling

These use the particle mass and temperature to determine velocity scale.

## Angular-momentum generation

The file includes angular-momentum distribution classes and helper functions used to initialise rotational state consistently with the supplied moments of inertia and temperature.

## Why this file matters

If a simulation starts with the wrong spatial structure, velocity scale, or rotational state, the issue is usually here or in the particle definitions passed in by `Initialiser.py`.

## When to edit this file

Edit `Distributions.py` if you want to:

- add a new position distribution
- add a new velocity distribution
- change how custom gradients are interpreted
- change how angular momentum is sampled