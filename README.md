# DED-Ecological-Modelling-code
This repository contains R code for a paper by Bajeux, Arino, Portet and Westwood, published in Ecological Modelling in 2020. The code here does not cover the entirety of the content of the paper. In particular, we do not provide specific code to generate all the figures in the paper, although they can easily be generated using the code here.

### Important
Some of the code requires substantial computing power (multiple cores) or RAM (>128 GB). This will be indicated in the function description below.

### R scripts
Although the `R` scripts are in the directory `CODE`, they are documented here. Most scripts have a prefixed name to indicate the stage at which they come into play.

- `pre_` are preprocessing scripts that typically need to be run before simulations can be performed.
- `run_` scripts run the simulations.
- `post_` scripts are run after simulations to process the results.
