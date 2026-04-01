# Assets and Folders

## `project/initialisation_files/`

This folder stores reusable text files describing simulation setups. These are the files shown when the programme asks you to choose a saved initialisation.

Examples already present in the repository include:

- water-only setups
- mixed-species setups
- atmosphere-style gradient setups

These files are produced by `Initialiser.py` when a manually entered configuration is saved.

## `project/plotted_graphs/`

This folder stores graph output. Each simulation run writes its PNG files into a subfolder named after the simulation.

Examples of subfolder names are usually based on:

- the initialisation filename, if a saved configuration was loaded
- a user-entered simulation name, if the run was created manually without saving

## `doc/project-file-guides/`

This folder contains the documentation added for the main source files so you can navigate the codebase by responsibility instead of by filename alone.