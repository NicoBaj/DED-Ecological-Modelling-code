# DED-Ecological-Modelling-code

### Simulation on data and neighbourhoods used in the article
To launch a simulation that uses the same dataset as the article (tree inventory dated 28 January 2020 and one of the two neighbourhoods PCP or NRH), run the script `run_one_sim.R` and set the logical gate `SIMULATIONS_ARTICLE` to `TRUE`. Then choose the neighbourhood, the type of initial conditions and the values of the main parameters.

### Simulation from newly updated tree inventory or with different neighbourhoods
To launch a new simulation, follow these steps: 
1. Run `pre_load_trees_and_compute_tree_heights.R` if you want to refresh the tree inventory. If you do, then you must also run `pre_roots_vs_routes.R` to eliminate connections between root systems intersected by roads. **Warning**: a minimum of **90 GB** of RAM is required to run the latter. The `DATA` directory contains files resulting from running these scripts, so you can try simulations even if you do not possess this amount of RAM (i.e., skip directly to step 2); they are named with the patterns `tree_inventory_elms_yyyy-mm-dd.Rds` and `elms_distances_roots_yyyy-mm-dd.Rds`, respectively.
3. Run `pre_neighbourhoods_proba_roots.R` with a selected neighbourhood.
4. Run `pre_network_beetles.R` with the chosen values for the beetle maximum dispersal distance `R_B` and the selected neighbourhood.
5. Run `run_one_sim.R` with `SIMULATIONS_ARTICLE = FALSE`, a value of `R_B` within the chosen values and the name of the neighbourhood.

### R scripts
All `R` scripts are documented here. Most scripts have a prefixed name to indicate the stage at which they come into play.

- `pre_` are preprocessing scripts that typically need to be run before simulations can be performed.
- `run_` scripts run the simulations.

#### Preprocessing 
- `set_directories.R` sets the directories.
- `pre_load_trees_and_compute_tree_heights.R` loads the tree inventory from the [City of Winnipeg Open Data portal](https://data.winnipeg.ca/), selects American Elms, computes the heights of each tree, extracts lat/lon information for each tree and saves the result for later use as file named `tree_inventory_elms_yyyy-mm-dd.Rds`, where `yyyy-mm-dd` is the date at which the query was performed.
- `pre_roots_vs_routes.R` (where route is intended as the French for road) loads Openstreetmap (OSM) data for the geography: roads, parking lots, railroads and rivers. It computes the distance between each pair of trees in the elms file then sets up line segments between all trees that are less than some threshold maximum distance away (we use 6 times the maximum tree height). It then removes from this graph all edges that intersect at least one of the OSM objects. The result is saved as `elms_distances_roots_yyyy-mm-dd.Rds`, where `yyyy-mm-dd` is the date of the elm tree inventory file used. Remark that we compute distances in the entire city. ***The amount of RAM required for some of the operations is therefore quite substantial. To run the code comfortably if other programs are running, 128 GB of RAM is recommended***. Note that you could adapt the code to operate on smaller areas, thereby reducing the need for RAM.
- `pre_neighbourhoods_proba_roots.R` loads the elm tree inventory created by `pre_load_trees_and_compute_tree_heights.R` and the roots-related dataset obtained in `pre_roots_vs_routes.R`. It selects trees from chosen neighbourhoods and saves this new tree inventory. Then, it selects, in the roots-related dataset, the rows for which the trees are in the neighbourhood and computes, for all pairs of trees in the neighbourhood, the probability to become infected through the root system. It saves this dataset.
- `pre_network_beetles.R` loads the tree inventory created by `pre_neighbourhoods_proba_roots.R` and saves, for the selected values of the beetle maximum dispersal distance `R_B`: the neighbours of each tree and information on neighbour trees.

#### Processing
- `functions_pre_simulations.R` has functions to set initial conditions and the simulation environment.
- `functions_simulation.R` contains functions required during one simulation.
- `run_one_sim.R` launches one simulation. Once the pre-processing is done, it selects the neighbourhood, the main parameter values and initial conditions. Then, it sources the functions required for the simulation.
