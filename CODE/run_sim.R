##run_one_sim.R: to run 

# File that conducts the main simulation run.
# Works by logical control through flow.

# LOGICAL GATES
#just a story of initial conditions ...
EVAL_IC = TRUE #do we need to create a new IC ?
EVAL_IC_RANDOM = FALSE #do we want random IC
PLOT_SIM = TRUE #if one sim is launched, do you want to see the DED spread evolution?

#choose the neighbourhood in which the simulation is run
Neighbourhood = "MIXED_PULBERRY_CRESCENT_PARK" #PCP
# Neighbourhood = "MIXED_PULBERRY_CRESCENT_PARK" #

#main parameters
input = list()
input$maxD = 20 #max distance that beetles fly during one time step
input$pr   = 0.5 #max proba for an infected tree to infect another one by root infection
input$pi   = 0.02 #proba that one beetle infects successfully one tree
input$sdt  = 0.98 #proba for beetles to survive one time step

#do we want to run the code with the root infection route ?
roots = TRUE

#number of simulations required for the set of parameters given above
nb_sims = 1

#if nb_sims>1, say if you want to run sims in parallel
RUN_PARALLEL = FALSE

# The three types of IC used in the paper
IC_type="cluster"
# IC_type = "2clusters"
# IC_type = "random"


IC_beetles = 500 #nb of inf beetles in each infected tree
IC_radius = 96 #radius of the cluster in m, to remove for 2clusters

if (IC_type == "2clusters"){ # this is for PCP
  IC_radius1 = 47.5
  IC_radius2 = 280
  IC_radius = list(r1=IC_radius1,r2=IC_radius2)
}
if(Neighbourhood=="MIXED_PULBERRY_CRESCENT_PARK" & IC_type == "cluster"){
  IC_radius = 96
}

# the next parameter can be an integer (classic approach) or a character chain that specifies the number of dead trees in function of the IC_radius for instance
IC_number_dead_trees = "ic_radius"
#old: "ic_radius"
# the following is used if random IC
if (EVAL_IC_RANDOM){
  if(Neighbourhood == "MIXED_PULBERRY_CRESCENT_PARK"){
    IC_number_dead_trees = 38
  }else if(Neighbourhood == "NORTH_RIVER_HEIGHTS"){
    IC_number_dead_trees = 53 ## to change that
  }
}

#Put here the initial and final dates for the simulation(s)
start_date = "2019-08-01"
end_date = "2020-12-31"

# if we need one particular initial condition, then put the file here
# to put the name of the file, go to sim_parallel...R and put the file name there
if (!EVAL_IC){
  dir_ic = "dir_file_IC"
}

# Libraries needed for all the processes
library(sqldf)
library(igraph)
library(Matrix)
library(R.utils)
library(markovchain)
library(parallel)
library(ISOweek)
library(stringr)
library(tictoc)
library(tgp)
library(poibin)

# Make code location aware using the library rprojroot
if(!require("rprojroot")){
  install.packages("rprojroot")
  require("rprojroot")
}
# Where is the current file? Needs to be sourced
this_file <- rprojroot::thisfile()
# Set code directory
TOP_DIR_CODE <- dirname(this_file)
# Set top directory
TOP_DIR <- dirname(TOP_DIR_CODE)

# Source code that sets directories. Always run. NO: The TOP_DIR_CODE
#source(sprintf("%s/set_directories.R", TOP_DIR_CODE))

# Source the pre processing (create env,loading files...)
source(sprintf("%s/pre_set_various.R",TOP_DIR_CODE))
# Source all the functions needed for the simulation
source(sprintf("%s/sim_functions.R",TOP_DIR_CODE))
# Source the functions that create initial condition
source(sprintf("%s/pre_IC.R",TOP_DIR_CODE))
# Source the R script for the simulation
source(sprintf("%s/sim.R",TOP_DIR_CODE))
