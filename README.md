# DED-Ecological-Modelling-code
This repository contains R code for a paper by Bajeux, Arino, Portet and Westwood, submitted to Ecological Modelling in 2020. The code here does not cover the entirety of the content of the paper. In particular, we do not provide specific code to generate all the figures in the paper, although they can easily be generated using the code here.

### Important
Some of the code requires substantial computing power (multiple cores) or RAM (>90 GB). This will be indicated in the function descriptions.

### Directories
-`CODE` contains all R scripts.

-`DATA` contains all data used for simulations (tree inventories, model parameter values, beetle <img src="https://render.githubusercontent.com/render/math?math=\mathcal{N}^B"> and root <img src="https://render.githubusercontent.com/render/math?math=\mathcal{N}^R"> networks).

-`RESULTS` will contain simulation results (data and figures).

-`simulation_DED_Shiny` is a barebones Shiny app allowing to run simple simulations. See that folder for details.

### Run simulations

#### Simulations on data and neighbourhoods used in the article
To launch one simulation using data and neighbourhoods from the article, go to `CODE/run_one_sim.R` and set `SIMULATIONS_ARTICLE` to `TRUE` and other inputs of the code (neighbourhood, initial conditions and main parameters). See `CODE/README.md` for more details.

#### Simulations on newly updated tree inventory 
To launch one simulation using newly updated data, see instructions in `CODE/README.md`.

#### Run simulations in the Shiny app
The file `simulation_DED_Shiny/app.R` allows to run extremely simple simulations in a Shiny app, i.e., in a easy to use graphical interface.
