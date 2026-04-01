# Initialiser.py

## Purpose

`Initialiser.py` is responsible for collecting all simulation input. It supports both interactive manual entry and loading from saved text files.

## Main responsibilities

- convert masses and moments of inertia between unit systems
- create and manage the initialisation-files folder
- save manual configurations to text files
- list and load saved configurations
- validate units and distribution-type selections
- build the particle-definition lists consumed by `Interface.py`

## Key functions

### `Initialise()`

This is the main entry point. It asks whether to use file-based initialisation.

- If the user chooses a file, it calls `ChooseFile()` and returns the parsed values.
- Otherwise it walks through manual input for the box, time settings, particle definitions, random seed, and optional save-to-file step.

### `ChooseFile()`

This shows the available files in `project/initialisation_files/`, asks the user which one to load, parses it, and returns the structured parameters.

### `SaveInputToFile(...)`

This writes a manual configuration back out as a reusable text file with comment headers.

### `InitialiseBox()`, `InitialiseTime()`, `InitialiseParticles()`

These functions handle the interactive prompts for each section of the configuration.

## Data shape returned for particles

Each particle type is represented as a list containing, in order:

1. particle type label
2. particle count
3. mass in kg
4. degrees of freedom
5. inertia and vibrational metadata dictionary
6. collision radius
7. initial temperature
8. position distribution type
9. velocity distribution type
10. optional gradient expression

That structure is later consumed by `Interface.py` and `Distributions.py`.

## Supported unit labels

Length units:

- `A`
- `nm`
- `mim`
- `mm`
- `cm`
- `m`

Time units:

- `ns`
- `mis`
- `ms`
- `s`

## Important behaviour

- Manual runs can be saved as reusable text files.
- File-based runs preserve the source filename so the simulation can inherit a matching output-folder name.
- Custom gradient position distributions can store a symbolic expression in the saved file.

## When to edit this file

Edit `Initialiser.py` if you want to change the input format, add new units, add new particle metadata fields, or change how initialisation files are saved and loaded.