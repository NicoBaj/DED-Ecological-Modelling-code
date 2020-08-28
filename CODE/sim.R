###############################################
## Spread of the Dutch Elm Disease
## Code to run simulations in parallel
###############################################

# Constants will go in a list called sim_constants.
sim_constants = list()

# Set part of the name for the result directory
sim_constants$sim_core = sim_core

# Set directories
sim_constants$DIRS = set_directories(TOP_DIR)

# Set the neighbourhood
sim_constants$Neighbourhood = Neighbourhood


# Read parameters from the csv file and replace those choosen in the main file (input)
sim_constants = read_parameters(sim_constants,input) # create default_params

# Root infection
sim_constants$roots = roots

# Set output directory
sim_constants = set_output_location(sim_constants)

# Set up simulations
sim_constants$time = set_sim_time(sim_constants)

# Nb of simulations
sim_constants$nb_sims = nb_sims

# Set up small demography matrices BY DEFAULT (they can change if we change parameters inside the matrices)
# sim_constants$matrices = set_demography_matrices(sim_constants$default_params) matrices are NOT constant then ...

# Set up other stuffs
sim_constants$other = set_other_constants(sim_constants)

#Set the gates
sim_constants$GATES = list()
sim_constants$GATES$PLOT_SIM = PLOT_SIM #add other gates if needed

# Parameters (and almost everything else) will go in a list called sim_params.
sims_params = list()

sims = list()
sims$params = set_sim_environment(sim_constants=sim_constants)

# Set the initial conditions
# sim_constants$initial_set_up = list()
# sim_constants$initial_set_up$IC_type=IC_type
# sim_constants$initial_set_up$IC_radius=IC_radius
# sim_constants$initial_set_up$IC_beetles=IC_beetles

if (IC_type=="random"){#if IC is random, then each sim has a different IC
  sims$IC = list()
  for (i in 1:nb_sims){
    sims$IC[[i]] = create_IC(sim_constants = sim_constants,Elms=sim_constants$default_params$elms,IC_type,IC_beetles,IC_radius,IC_number_dead_trees)
  }
}else{#if IC is cluster or 2clusters
  sims$IC = create_IC(sim_constants = sim_constants,Elms=sim_constants$default_params$elms,IC_type,IC_beetles,IC_radius,IC_number_dead_trees=0) #put IC_number_dead_trees to 0 here, the function will replace its value by the good nb of initially infected trees
}
sim_constants$default_params$IC= sims$IC

# sim_constants$initial_set_up$IC_number_dead_trees=sims$IC$IC_number_dead_trees

if (FALSE){#if needed, we can create an IC using the functions in pre_IC.R and read the file directly
  sims$IC = readRDS("name_file_IC.rds")
}


for (i in 1:sim_constants$nb_sims) {
  sims_params[[i]] = list()
  sims_params[[i]]$params = sims$params$params[[i]]
  if(IC_type=="random"){ #for random IC: there is a different IC for each sim
    sims_params[[i]]$IC = sims$IC[[i]] 
  }else{ #for cluster and 2clusters IC: there is only one IC
    sims_params[[i]]$IC = sims$IC 
  }
  
}

# Run simulations
if (RUN_PARALLEL) {
  # Detect number of cores, use all but 1
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
  results = parLapply(cl = cl, X = sims_params, fun =  function(x) system_over_time(x,sim_constants))
  # Stop cluster
  stopCluster(cl)
  timeProcessing=tictoc::toc()
} else {
  tictoc::tic()
  # Run computation
  results = lapply(X = sims_params, FUN = function(x) system_over_time(x,sim_constants))
  timeProcessing=tictoc::toc()
}




#########################################################################
## Here, we save all the results in one file
#########################################################################


## We just need to save the output results: tree states and beetle populations in results[[i]]$sim_output
if (SIM_SAVE){
  tree_states = list()
  beetles = list()
  paramet = list()
  vec_year = seq(1,length(sim_constants$time$idx),by=53)
  
  for (i in 1:sim_constants$nb_sims){
    tree_states[[i]] = results[[i]]$sim_output$status_trees[,vec_year]
    beetles[[i]]     = results[[i]]$sim_output$matPopByTrees
    paramet[[i]]     = results[[i]]$sim_param$params
  }
  
  saveRDS(sims$IC, 
          file = sprintf("%s/ic.Rds",sim_constants$output_dir))
  
  saveRDS(tree_states, 
          file = sprintf("%s/tree_states_maxD%03g_pb_inf%03g_pr%03g.Rds",
                         sim_constants$output_dir,
                         sim_constants$default_params$maxD,
                         sim_constants$default_params$proba_infection*100,
                         sim_constants$default_params$p_r*100))
  
  saveRDS(beetles, 
          file = sprintf("%s/beetles_maxD%03g_pb_inf%03g_pr%03g.Rds",
                         sim_constants$output_dir,
                         sim_constants$default_params$maxD,
                         sim_constants$default_params$proba_infection*100,
                         sim_constants$default_params$p_r*100))  
  
  # If needed, save the parameters as well (can be a large file)
  if(FALSE){
    saveRDS(paramet,
            file = sprintf("%s/parameters_maxD%03g_pb_inf%03g_pr%03g.Rds",
                           sim_constants$output_dir,
                           sim_constants$default_params$maxD,
                           sim_constants$default_params$proba_infection*100,
                           sim_constants$default_params$p_r*100))
  }
  
  
}
