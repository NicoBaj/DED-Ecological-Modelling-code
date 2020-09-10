##run_one_sim.R: 
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
SIM_SAVE = FALSE #do we save the outputs ?

#inputs are obtained from OUT
OUT = RESULTS_NEW()
input = OUT # Main model parameter inputs are defined here

# 1- Neighbourhood
if(input$NGH == "North River Heights (NRH)"){
  Neighbourhood = "NORTH_RIVER_HEIGHTS"
}else if(input$NGH == "Mixed Pulberry-Crescent Park (PVC)"){
  Neighbourhood = "MIXED_PULBERRY_CRESCENT_PARK"
}
date_TI_file = "2020-01-28"

#2- IC
if(input$IC=="Cluster"){
  IC_type="cluster"
}else if(input$IC=="Two clusters"){
  IC_type = "2clusters"
}else if(input$IC=="Random"){
  IC_type = "random"
}

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
end_date = "2030-12-31"

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

# Set directories
source("set_directories.R")
# Set save directory to include date of data file
DIRS$nbhd_and_date = sprintf("%s%s", DIRS$prefix_data_date, date_TI_file)
# Set directory for saving in this script
DIRS$preproc_dists = sprintf("%s/%s", DIRS$nbhd_and_date, DIRS$suffix_preproc_dists)


# Source the functions that help set up the simulation
source(sprintf("%s/functions_pre_simulation.R", DIRS$CODE))
# Source all the functions needed for the simulation
source(sprintf("%s/functions_simulation.R", DIRS$CODE))


## MAIN SIMULATION PART
# Constants will go in a list called sim_constants.
sim_constants = list()
#Set the gates
sim_constants$GATES = list()
sim_constants$GATES$PLOT_SIM = PLOT_SIM #add other gates if needed
sim_constants$GATES$SIMULATIONS_ARTICLE = SIMULATIONS_ARTICLE

# Add directories to constants for easy use
sim_constants$DIRS = DIRS

# Set the neighbourhood
sim_constants$neighbourhood = neighbourhood

# Set part of the name for the result directory
if (SIMULATIONS_ARTICLE) {
  sim_constants$sim_core = sim_core
} else {
  nb_trees = readRDS(sprintf("%s/sim_core_%s.Rds",
                             sim_constants$DIRS$nbhd_and_date,
                             neighbourhood))
  sim_constants$sim_core = sprintf("%strees",nb_trees)
}

# Read parameters from the csv file and replace those choosen in the main file (input)
sim_constants = read_parameters(sim_constants, input) # create default_params
# Root infection
sim_constants$roots = roots
# Set output directory
sim_constants = set_output_location(sim_constants)
# Set up simulations
sim_constants$time = set_sim_time(sim_constants)
# Nb of simulations
sim_constants$nb_sims = nb_sims

# Set up small demography matrices BY DEFAULT (they can change if we change parameters inside the 
# matrices)
# sim_constants$matrices = set_demography_matrices(sim_constants$default_params) matrices are NOT constant then ...

# Set up other stuffs
sim_constants$other = set_other_constants(sim_constants)

# Parameters (and almost everything else) will go in a list called sim_params.
sims_params = list()

sims = list()
sims$params = set_sim_environment(sim_constants=sim_constants)

# Set the initial conditions
# sim_constants$initial_set_up = list()
# sim_constants$initial_set_up$IC_type=IC_type
# sim_constants$initial_set_up$IC_radius=IC_radius
# sim_constants$initial_set_up$IC_beetles=IC_beetles

if (IC_type=="random") {
  #if IC is random, then each sim has a different IC
  sims$IC = list()
  for (i in 1:nb_sims) {
    sims$IC[[i]] = create_IC(sim_constants = sim_constants,
                             elms = sim_constants$default_params$elms,
                             IC_type,
                             IC_beetles,
                             IC_radius,
                             IC_number_dead_trees)
  }
} else {
  #if IC is cluster or 2clusters
  sims$IC = create_IC(sim_constants = sim_constants,
                      elms = sim_constants$default_params$elms,
                      IC_type,
                      IC_beetles,
                      IC_radius,
                      IC_number_dead_trees = 0) #put IC_number_dead_trees to 0 here, the function will replace its value by the good nb of initially infected trees
}
sim_constants$default_params$IC= sims$IC

# sim_constants$initial_set_up$IC_number_dead_trees=sims$IC$IC_number_dead_trees

if (FALSE) {
  #if needed, we can create an IC using the functions in pre_IC.R and read the file directly
  sims$IC = readRDS("name_file_IC.rds")
}


for (i in 1:sim_constants$nb_sims) {
  sims_params[[i]] = list()
  sims_params[[i]]$params = sims$params$params[[i]]
  if(IC_type=="random"){ #for random IC: there is a different IC for each sim
    sims_params[[i]]$IC = sims$IC[[i]] 
  } else { #for cluster and 2clusters IC: there is only one IC
    sims_params[[i]]$IC = sims$IC 
  }
}

# Run simulations
if (RUN_PARALLEL) {
  # Detect number of cores
  no_cores <- detectCores()
  # Initiate cluster
  tictoc::tic()
  cl <- makeCluster(no_cores)
  # Export needed variables
  clusterEvalQ(cl,{
    library(markovchain)
    library(Matrix)
    library(poibin)
  })
  clusterExport(cl,
                c("system.over.time",
                  "convert_state_to_number",
                  "merge_updates",
                  "mvt.beetles",
                  "proba.distance",
                  "demography.matrices",
                  "sim_constants",
                  "new.transition.matrices"),
                envir = .GlobalEnv)
  # Run computation
  results = parLapply(cl = cl, X = sims_params, 
                      fun = function(x) system_over_time(x, sim_constants))
  # Stop cluster
  stopCluster(cl)
  timeProcessing=tictoc::toc()
} else {
  # Run computation
  results = lapply(X = sims_params, FUN = function(x) system_over_time(x, sim_constants))
}




#########################################################################
## Save all the results in one file
#########################################################################


## We just need to save the output results: tree states and beetle populations in results[[i]]$sim_output
if (SIM_SAVE) {
  tree_states = list()
  beetles = list()
  paramet = list()
  vec_year = seq(from = 1, to = length(sim_constants$time$idx), by = 53)
  for (i in 1:sim_constants$nb_sims) {
    tree_states[[i]] = results[[i]]$sim_output$status_trees[,vec_year]
    beetles[[i]]     = results[[i]]$sim_output$matPopByTrees
    paramet[[i]]     = results[[i]]$sim_param$params
  }
  saveRDS(sims$IC, 
          file = sprintf("%s/ic.Rds",sim_constants$DIRS$output_dir))
  saveRDS(tree_states, 
          file = sprintf("%s/tree_states_Rb%03g_pb_inf%03g_pr%03g.Rds",
                         sim_constants$DIRS$output_dir,
                         sim_constants$default_params$R_B,
                         sim_constants$default_params$p_i*100,
                         sim_constants$default_params$p_r*100))
  saveRDS(beetles, 
          file = sprintf("%s/beetles_Rb%03g_pb_inf%03g_pr%03g.Rds",
                         sim_constants$DIRS$output_dir,
                         sim_constants$default_params$R_B,
                         sim_constants$default_params$p_i*100,
                         sim_constants$default_params$p_r*100))  
  
  # If needed, save the parameters as well (can be a large file)
  if (FALSE) {
    saveRDS(paramet,
            file = sprintf("%s/parameters_Rb%03g_pb_inf%03g_pr%03g.Rds",
                           sim_constants$DIRS$output_dir,
                           sim_constants$default_params$R_B,
                           sim_constants$default_params$p_i*100,
                           sim_constants$default_params$p_r*100))
  }
}
