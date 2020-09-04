# DED-Ecological-Modelling-code
This repository contains R code for a paper by Bajeux, Arino, Portet and Westwood, published in Ecological Modelling in 2020. The code here does not cover the entirety of the content of the paper. In particular, we do not provide specific code to generate all the figures in the paper, although they can easily be generated using the code here.

### Important
Some of the code requires substantial computing power (multiple cores) or RAM (>128 GB). This will be indicated in the function description below.

### R scripts
Although the `R` scripts are in the directory `CODE`, they are documented here. Most scripts have a prefixed name to indicate the stage at which they come into play.

- `pre_` are preprocessing scripts that typically need to be run before simulations can be performed.
- `run_` scripts run the simulations.
- `post_` scripts are run after simulations to process the results.

#### Preprocessing 
- `pre_load_trees_and_compute_tree_heights.R` loads the tree inventory from the [City of Winnipeg Open Data Portal](), selects American Elms, computes the heights of each tree, extracts lat/lon information for each tree and saves the result for later use.
- `pre_roots_vs_routes.R` (where route is intended as the French for road) loads Openstreetmap (OSM) data for the geography: roads, parking lots, railroads and rivers. It computes the distance between each pair of trees in the elms file then sets up line segments between all trees that are less than some threshold maximum distance away (we use 6 times the tree heights). It then removes from this graph all edges that intersect at least one of the OSM objects. Remark that we compute distances in the entire city. The amount of RAM required for some of the operations is therefore quite substantial. To run the code comfortably if other programs are running, 128 GB of RAM is therefore recommended.
- `pre_neighbourhoods_proba_roots.R` loads the elm tree inventory created in `pre_load_trees_and_compute_tree_heights.R` and the roots-related dataset obtained in `pre_roots_vs_routes.R`. It selects trees from chosen neighbourhoods and saves this new tree inventory. Then, it selects in the roots-related dataset the rows for which the trees are in the neighbourhood. It then computes, for all couples of trees in the neighbourhood, the probability to become infected through the root system. It saves this dataset as "Proba_roots_Neighbourhood.Rds".
-`pre_network_beetles.R` loads one tree inventory of a given neighbourhood and returns, for the selected values of the beetle maximum dispersal distance $R_B$, three things (all are lists): 1) the neighbours of 
each tree, 2) their position in the neighbourhood and 3) their distance.


-`run_sim.R` launches one simulation. Once the pre-processing done, it selects the neighbourhood, the main parameter values and initial conditions. Then, it sources the functions required for the simulation.
- `pre_functions_IC.R` stores functions to set initial conditions.
- `pre_set_various.R` sets the environment to run simulations.
 
 
- `sim_functions.R` stores functions required during one simulation.
- `sim.R` is the code that runs 
