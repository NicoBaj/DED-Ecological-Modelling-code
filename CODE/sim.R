###############################################
## Spread of the Dutch Elm Disease
## Code to run simulations in parallel
###############################################

# Constants will go in a list called sim_constants.
sim_constants = list()

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

# Read parameters
sim_constants = read_parameters(sim_constants,input) # create default_params
# Root infection
sim_constants$roots = TRUE
# Set output directory
sim_constants = set_output_location(sim_constants)
# Set up simulations
sim_constants$time = set_sim_time(sim_constants)

# Nb of simulations
sim_constants$nb_sims = 1

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

# Parameters (and almost everything else) will go in a list called sim_params.
sims_params = list()

sims = list()
sims$params = set_sim_environment(sim_constants=sim_constants)

# Choice of initial conditions
# new_IC = TRUE
if(EVAL_IC){
  sims$IC = create_IC(sim_constants = sim_constants,Elms=sim_constants$default_params$elms,IC_type,IC_beetles,IC_radius,IC_number_dead_trees)
  dir_ic=sprintf("%s/PARAMS/",TOP_DIR_DATA_OUTPUT)
  if(is.list(IC_radius)){
    saveRDS(sims$IC,file = sprintf("%sIC_%s_radiuses_%s_%s.RData",dir_ic,sim_constants$Neighbourhood,IC_radius$r1,IC_radius$r2))
  }else{
    saveRDS(sims$IC,file = sprintf("%sIC_%s_radius_%s.RData",dir_ic,sim_constants$Neighbourhood,IC_radius))
  }
  
  sim_constants$default_params$IC= sims$IC
}else if (EVAL_IC_RANDOM){
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

sims_params[[1]] = list()
sims_params[[1]]$params = sims$params$params[[1]]
sims_params[[1]]$IC = sims$IC # only initial beetles, initial tree status and initial demography

# Run simulations
results = lapply(X = sims_params, FUN = function(x) system.over.time(x,sim_constants))


# saveRDS(sims$IC, file = sprintf("%s/ic.RData",sim_constants$output_dir))

#########################################################################
## Here, we save all the results in one file
#########################################################################


## We just need to save the output results (TREE STATES here) which is in results[[i]]$sim_output
tree_states = list()
beetles = list()
paramet = list()
root_or_beetle = list()
vec_year = seq(1,length(sim_constants$time$idx),by=53)
  tree_states[[1]] = results[[1]]$sim_output$status_trees[,vec_year]
  beetles[[1]]     = results[[1]]$sim_output$matPopByTrees
  paramet[[1]]     = results[[1]]$sim_param$params
  root_or_beetle[[1]] = results[[1]]$root_or_vec

saveRDS(sims$IC, file = sprintf("%s/ic.RData",sim_constants$output_dir))
saveRDS(tree_states, 
        file = sprintf("%s/tree_states_maxD%03g_pb_inf%03g_pr%03g.RData",
                       sim_constants$output_dir,
                       sim_constants$default_params$maxD,
                       sim_constants$default_params$proba_infection*100,
                       sim_constants$default_params$p_r*100))
# saveRDS(beetles, file = sprintf("%s/beetles%05s_maxD%03g_pb_inf%03g_pr%03g.RData",sim_constants$output_dir,sim_constants[["fileNb"]],sim_constants$default_params$maxD,sim_constants$default_params$proba_infection*100,sim_constants$default_params$p_r*100))  
saveRDS(root_or_beetle, 
        file = sprintf("%s/root_or_beetle_maxD%03g_pb_inf%03g_pr%03g.RData",
                       sim_constants$output_dir,
                       sim_constants$default_params$maxD,
                       sim_constants$default_params$proba_infection*100,
                       sim_constants$default_params$p_r*100))
# saveRDS(sim_constants$default_params$maxD, file = sprintf("%s/maxD",sim_constants$output_dir))

# remove the elms and the preproc from sim_constants:
# These can be accessed via sim_constants$FILES[[i]]
sim_constants$default_params$elms<-NULL
sim_constants$default_params$preproc <- NULL

states = results[[1]]$sim_output$status_trees
states = convert_state_to_number(states)
elms = readRDS(sprintf("%s/Elms_Neighbourhood/Elms_MIXED_PULBERRY_CRESCENT_PARK.RData",TOP_DIR_DATA))
