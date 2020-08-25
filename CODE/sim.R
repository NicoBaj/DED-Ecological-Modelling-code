###############################################
## Spread of the Dutch Elm Disease
## Code to run simulations in parallel
###############################################

# Constants will go in a list called sim_constants.
sim_constants = list()

#################
##NOTE_NB: I need to remove that once we change the way we save the pre_processing stuff
#################
# Set part of the name for the result directory
sim_constants$sim_core = "1513Trees"

# Set directories
sim_constants$DIRS = set_directories(TOP_DIR)

# Set the neighbourhood
sim_constants$Neighbourhood = Neighbourhood
# Set the important initial parameters
sim_constants$initial_set_up = list()
sim_constants$initial_set_up$IC_type=IC_type
sim_constants$initial_set_up$IC_radius=IC_radius
sim_constants$initial_set_up$IC_beetles=IC_beetles
sim_constants$initial_set_up$IC_number_dead_trees=IC_number_dead_trees

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

#if we want to compare random IC
if(EVAL_IC_RANDOM){
  if(is.null(sim_constants$fileNb)){
    sim_constants$fileNb = 1
  }
  # sim_constants$default_params$maxD = set_maxD(sim_constants,sim_constants$fileNb)# ten values
  sim_constants$default_params$maxD = set_maxD_IC_random(sim_constants,sim_constants$fileNb) #only 3 values
}

#Set the gates
sim_constants$GATES = list()
sim_constants$GATES$PLOT_SIM = PLOT_SIM #add other gates if needed

# Parameters (and almost everything else) will go in a list called sim_params.
sims_params = list()

sims = list()
sims$params = set_sim_environment(sim_constants=sim_constants)

# Choice of initial conditions
# 
if(EVAL_IC){
  sims$IC = create_IC(sim_constants = sim_constants,Elms=sim_constants$default_params$elms,IC_type,IC_beetles,IC_radius,IC_number_dead_trees)
  if(is.list(IC_radius)){#in the case of two clusters
    saveRDS(sims$IC,file = sprintf("%s/IC_%s_radiuses_%s_%s.RData",sim_constants$DIRS$DATA,sim_constants$Neighbourhood,IC_radius$r1,IC_radius$r2))
  }else{
    saveRDS(sims$IC,file = sprintf("%s/IC_%s_radius_%s.RData",sim_constants$DIRS$DATA,sim_constants$Neighbourhood,IC_radius))
  }
  sim_constants$default_params$IC= sims$IC
}else if (EVAL_IC_RANDOM){#if IC is random, then each sim has a different IC
  sims$IC = list()
  for (i in 1:nb_sims){
    sims$IC[[i]] = create_IC(sim_constants = sim_constants,Elms=sim_constants$default_params$elms,IC_type,IC_beetles,IC_radius,IC_number_dead_trees)
    dir_ic=sprintf("%s/PARAMS/",TOP_DIR_DATA_OUTPUT)
    sim_constants$default_params$IC= sims$IC
  }
  saveRDS(sims$IC,file = sprintf("%sIC_%s_random.RData",dir_ic,sim_constants$Neighbourhood))
}else{
  sims$IC = readRDS(sprintf("%sIC_%s_radius_%s.RData",dir_ic,sim_constants$Neighbourhood,IC_radius))
}


for (i in 1:sim_constants$nb_sims) {
  sims_params[[i]] = list()
  sims_params[[i]]$params = sims$params$params[[i]]
  sims_params[[i]]$IC = sims$IC # only initial beetles, initial tree status and initial demography
}

# sims_params[[1]]        = list()
# sims_params[[1]]$params = sims$params$params[[1]]
# sims_params[[1]]$IC     = sims$IC # only initial beetles, initial tree status and initial demography

# Run simulations and plot in the same time the spread of the disease
# results = lapply(X = sims_params, FUN = function(x) system.over.time(x,sim_constants))


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
  results = parLapply(cl = cl, X = sims_params, fun =  function(x) system.over.time(x,sim_constants))
  # Stop cluster
  stopCluster(cl)
  timeProcessing=tictoc::toc()
} else {
  tictoc::tic()
  # Run computation
  results = lapply(X = sims_params, FUN = function(x) system.over.time(x,sim_constants))
  timeProcessing=tictoc::toc()
}




#########################################################################
## Here, we save all the results in one file
#########################################################################


## We just need to save the output results: tree states and beetle populations in results[[i]]$sim_output
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
        file = sprintf("%s/ic.RData",sim_constants$output_dir))

saveRDS(tree_states, 
        file = sprintf("%s/tree_states_maxD%03g_pb_inf%03g_pr%03g.RData",
                       sim_constants$output_dir,
                       sim_constants$default_params$maxD,
                       sim_constants$default_params$proba_infection*100,
                       sim_constants$default_params$p_r*100))

saveRDS(beetles, 
        file = sprintf("%s/beetles_maxD%03g_pb_inf%03g_pr%03g.RData",
                       sim_constants$output_dir,
                       sim_constants$default_params$maxD,
                       sim_constants$default_params$proba_infection*100,
                       sim_constants$default_params$p_r*100))  

# If needed, save the parameters as well (can be a large file)
if(FALSE){
  saveRDS(paramet,
          file = sprintf("%s/parameters_maxD%03g_pb_inf%03g_pr%03g.RData",
                         sim_constants$output_dir,
                         sim_constants$default_params$maxD,
                         sim_constants$default_params$proba_infection*100,
                         sim_constants$default_params$p_r*100))
}

