##run_sim.R: 
# File that conducts the main simulation run.

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


# LOGICAL GATES
SIMULATIONS_ARTICLE = TRUE # if true, launch simulations with tree inventory used in the article 
                           # (dated from 28th January 2020), else launch new simulation 
                           # (requires pre-processings)
PLOT_SIM = TRUE #if sims are launched, do you want to see the DED spread evolution?
SIM_SAVE = TRUE #do we save the outputs ?

#Inputs to define the type of simulations:
#1- choose the neighbourhood
#2- choose the initial condition (IC)
#3- choose the values of the main parameters

# 1- neighbourhood
# Neighbourhoods from the articles are "MIXED_PULBERRY_CRESCENT_PARK" (PCP) and "NORTH_RIVER_HEIGHTS" (NRH)
if (SIMULATIONS_ARTICLE) { # (un)comment to choose a neighbourhood
  # neighbourhood = "MIXED_PULBERRY_CRESCENT_PARK" #PCP
  neighbourhood = "NORTH_RIVER_HEIGHTS" #NRH
  date_TI_file = "2020-01-28"
} else {
  neighbourhood = "your_pick" # and replace spaces by "_" in the neighbourhood name
  date_TI_file = "2020-06-28" # Or something else..
}


#2- IC: uncomment your choice
IC_type = "cluster"
# IC_type = "2clusters"
# IC_type = "random"

#3- Main model parameters: update the values 
input = list()
input$R_B   = 100 #max distance that beetles fly during one time step. If using data from article, then choose from 20 to 380 by steps of 40. 
input$p_r   = 0.5 #max proba for an infected tree to infect another one by root infection
input$p_i   = 0.02 #proba that one beetle infects successfully one tree
input$s_dt  = 0.98 #proba for beetles to survive one time step

#
#do we want to run the code with the root infection route ?
roots = TRUE

#number of simulations required for the set of parameters given above
nb_sims = 1

#if nb_sims>1, say if you want to run sims in parallel
RUN_PARALLEL = FALSE

IC_beetles = 500 #nb of inf beetles in each infected tree

#Put here the initial and final dates for the simulation(s)
start_date = "2019-08-01"
end_date = "2021-12-31"

###############################################################
#Set up the simulation in function of choices in 1- and 2-:
if(SIMULATIONS_ARTICLE){
  if (neighbourhood == "MIXED_PULBERRY_CRESCENT_PARK") {
    sim_core = "1513Trees"
    if (IC_type == "cluster") {
      IC_radius = 96
    } else if (IC_type == "2clusters") {
      IC_radius1 = 47.5
      IC_radius2 = 280
      IC_radius = list(r1 = IC_radius1, r2 = IC_radius2)
      #These values of radii are done to get the same nb of infected trees at the initial time than in one cluster (96m)
    } else {
      IC_radius = 96
      IC_number_dead_trees = 38
    }
  } else if (neighbourhood == "NORTH_RIVER_HEIGHTS") {
    sim_core = "2004Trees"
    if (IC_type == "cluster") {
      IC_radius = 100
    } else if (IC_type == "2clusters") {
      IC_radius1 = 84
      IC_radius2 = 84
      IC_radius = list(r1 = IC_radius1, r2 = IC_radius2)
      #These values of radii are done to get the same nb of infected trees at the initial time than in one cluster (100m)
    }else{
      IC_radius = 100
      IC_number_dead_trees = 50
    }
  }
} else {#if this is a new simulation, put whatever you want in the following
  if (IC_type == "cluster") {
    IC_radius = 100
  } else if (IC_type == "2clusters") {
    IC_radius1 = 84
    IC_radius2 = 84
    IC_radius = list(r1=IC_radius1, r2=IC_radius2)
    #These values of radii are done to get the same nb of infected trees at the initial time than in one cluster (100m)
  } else {
    IC_radius = 100
    IC_number_dead_trees = 50
  }
}

# Set directories
source(sprintf("%s/CODE/set_directories.R", here::here()))

# Source the pre processing (create env,loading files...)
source(sprintf("%s/set_various.R", DIRS$CODE))
# Source all the functions needed for the simulation
source(sprintf("%s/functions_simulation.R", DIRS$CODE))
# Source the functions that create initial condition
source(sprintf("%s/functions_IC.R", DIRS$CODE))
# Source the R script for the simulation
source(sprintf("%s/sim.R", DIRS$CODE))
