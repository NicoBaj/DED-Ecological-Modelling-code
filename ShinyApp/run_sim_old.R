##run_one_sim.R: 
# File that conducts the main simulation run.

# print(commandArgs())
# 
# # Get output directory, simulation number and number of simulations in a file
# args <- commandArgs(trailingOnly = TRUE)
# # # Incorporate arguments read from the command line
# if (length(args)>0){
#   list.args=list()
#   for (i in 1:length(args)) {
#     strs = strsplit(args[i],"=")
#     list.args[[strs[[1]][1]]] = strs[[1]][2]
#   }
# }
# print(list.args)

# LOGICAL GATES
# PLOT_SIM = TRUE #if sims are launched, do you want to see the DED spread evolution?
SIM_SAVE = FALSE #do we save the outputs ?

#The three following inputs are crucial for the type of simulations:
#1- choose the neighbourhood
#2- choose the initial condition (IC)
#3- choose the values of the main parameters
#inputs are obtained from OUT
OUT = RESULTS_NEW()
input = OUT

PLOT_SIM=TRUE
PLOT_PROP=FALSE

# 1- Neighbourhood
if(input$NGH == "North River Heights (NRH)"){
  Neighbourhood = "NORTH_RIVER_HEIGHTS"
}else if(input$NGH == "Mixed Pulberry-Victoria Crescent-Crescent Park (PVC)"){
  Neighbourhood = "MIXED_PULBERRY_CRESCENT_PARK"
}
#add ngh if more ngh in the selectInput


#2- IC
if(input$IC=="Cluster"){
  IC_type="cluster"
}else if(input$IC=="Two clusters"){
  IC_type = "2clusters"
}else if(input$IC=="Random"){
  IC_type = "random"
}

#3- main parameters

#Set up the simulation in function of choices in 2- and 3-:
if(Neighbourhood=="MIXED_PULBERRY_CRESCENT_PARK"){
  sim_core = "1513Trees"
  if(IC_type == "cluster"){
    IC_radius = 96
  }else if(IC_type == "2clusters"){
    IC_radius1 = 47.5
    IC_radius2 = 280
    IC_radius = list(r1=IC_radius1,r2=IC_radius2)
    #These values of radii are done to get the same nb of infected trees at the initial time than in one cluster (96m)
  }else{
    IC_radius = 96
    IC_number_dead_trees = 38
  }
}else if (Neighbourhood == "NORTH_RIVER_HEIGHTS"){
  sim_core = "2004Trees"
  if(IC_type == "cluster"){
    IC_radius = 100
  }else if(IC_type == "2clusters"){
    IC_radius1 = 84
    IC_radius2 = 84
    IC_radius = list(r1=IC_radius1,r2=IC_radius2)
    #These values of radii are done to get the same nb of infected trees at the initial time than in one cluster (100m)
  }else{
    IC_radius = 100
    IC_number_dead_trees = 50
  }
}

#do we want to run the code with the root infection route ?
roots = TRUE

#number of simulations required for the set of parameters given above
nb_sims = 1

#if nb_sims>1, say if you want to run sims in parallel
RUN_PARALLEL = FALSE

IC_beetles = 500 #nb of inf beetles in each infected tree

#Put here the initial and final dates for the simulation(s)
start_date = "2019-08-01"
end_date = "2030-12-31"

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
source(sprintf("%s/pre_set_various.R",TOP_DIR_CODE),local = TRUE)
# Source all the functions needed for the simulation
source(sprintf("%s/sim_functions.R",TOP_DIR_CODE),local = TRUE)
# Source the functions that create initial condition
source(sprintf("%s/pre_IC.R",TOP_DIR_CODE),local = TRUE)
# Source the R script for the simulation
source(sprintf("%s/sim.R",TOP_DIR_CODE),local = TRUE)
