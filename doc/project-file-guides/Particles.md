# Particles.py

## Purpose

`Particles.py` contains both the per-particle data model and the vectorised array engine used to update many particles efficiently.

## `particle` class

The `particle` class stores the state of a single particle, including:

- index
- mass
- degrees of freedom
- moments of inertia
- vibrational modes
- collision radius
- position
- velocity
- angular momentum
- collision and impulse counters

It also provides scalar methods for:

- updating position
- checking collisions in the non-vectorised style
- computing kinetic energy
- computing rotational kinetic energy
- computing vibrational energy
- estimating temperature and pressure contributions

## `ParticleArray` class

This is the performance-critical part of the file. It mirrors the particle list into NumPy arrays so the simulation can work on many particles at once.

### What it stores

- masses, radii, particle types, and degrees of freedom
- positions, velocities, and angular momenta
- moments of inertia and vibrational-mode metadata
- accumulated collision and impulse counters

### What it does

- updates positions vectorially
- resolves particle-particle collisions using a spatial cell list
- handles wall collisions in each axis
- computes kinetic, rotational, and vibrational energies
- computes temperatures and pressure estimates
- computes pairwise distances in blocks to avoid forming a full dense `N x N` matrix
- computes nearest-neighbour distances
- syncs array state back to the underlying particle objects after the run

## Why this file matters

This file controls most of the runtime cost and memory behaviour of the simulation. Changes here can strongly affect correctness and performance.

## When to edit this file

Edit `Particles.py` when working on collision handling, performance, vectorisation, particle properties, or any calculation that needs direct access to all particle states.