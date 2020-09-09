# DATA directory

This directory contains data used in the simulations. Tree inventories, roots-related pre-processed networks and parameters are stored in this folder. Beetles networks are stored in the folder `/Preprocessing`.

#### Data used in the article

- Beetle networks are stored in the folder `/Preprocessing/article_pre_processing/`.

- Tree inventories from 28th January 2020 are in provided directly in `/DATA` with the date at the end of their name.

- Model parameters are in `parameters.csv`.

- Root networks are in the files `Proba_roots_*.Rds` in `/DATA`.

#### Data used for a new simulation

- After running the pre-processing scripts `pre_load_trees_and_compute_tree_heights.R` and `pre_roots_vs_routes.R`, it stores a new tree inventory and a dataset `elms_distances_roots*.Rds` in `/DATA`.

- After running `pre_neighbourhoods_proba_roots.R`, a dataset `Proba_roots_*.Rds` is stored in `/DATA/Roots`.

- After running `pre_network_beetles.R`, beetles-related networs are stored in `/DATA/Preprocessing/new_pre_processing/`.

- Model parameters are in `parameters.csv`.

Outputs of the simulation are stored in a new folder; three files are created:

- `beetles_*.Rds` stores the beetle population sizes in all trees at each week of the simulation.

- `tree_states_*.Rds` stores the status of trees at each year of the simulation.

- `ic.Rds`is the initial condition used for the simulation.


We did not package such results in the repository, so this file is here only to ensure that the directory is generated when you clone the repository.

