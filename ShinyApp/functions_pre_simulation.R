##pre_IC.R:
#functions that permit to create the three types of initial conditions (cluster, 2clusters and random)

###FIND_CENTER
#
#center of the map
find_center = function(elms){
  min_X = min(elms$X)
  max_X = max(elms$X)
  min_Y = min(elms$Y)
  max_Y = max(elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  return(center)
}

###SPECIAL_CENTER
#
#set the center of the third quadrant (bottom left)
special_center = function(elms){
  min_X = min(elms$X)
  max_X = max(elms$X)
  min_Y = min(elms$Y)
  max_Y = max(elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  third_quarter_center = list(x=(center$x+min_X)/2,y=(center$y+min_Y)/2)
  return(third_quarter_center)
}

###SPECIAL_CENTER2
#
#set the right center of the map (middle right)
special_center2 = function(elms){
  min_X = min(elms$X)
  max_X = max(elms$X)
  min_Y = min(elms$Y)
  max_Y = max(elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  middle_right_center = list(x=(center$x+max_X)/2,y=(max_Y+min_Y)/2)
  return(middle_right_center)
}

###SPECIAL_CENTER3
#
#set the center of the first quadrant (top right)
special_center3 = function(elms){
  min_X = min(elms$X)
  max_X = max(elms$X)
  min_Y = min(elms$Y)
  max_Y = max(elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  first_quadrant_center = list(x=(center$x+max_X)/2,y=(center$y+max_Y)/2)
  return(first_quadrant_center)
}

convert_state_to_number = function(M){
  res = mat.or.vec(dim(M)[1],dim(M)[2])
  res[which(M=="H")] = 1
  res[which(M=="S_W")] = 2
  res[which(M=="I_W")] = 3
  res[which(M=="S_D")] = 4
  res[which(M=="I_D")] = 5
  return(res)
}

###INITIALLY_INFECTED_TREES
#
#Set up the stages of trees at the initial time
initially_inf_trees = function(default_params, dead_sampling){
  current_stages = rep(1,default_params$N) #1 is for healthy
  initially_dead_trees_idx = dead_sampling 
  current_stages[initially_dead_trees_idx]=5 #5 is for dead and infected tree
  
  current_stages[which(current_stages==1)] = "H"
  current_stages[which(current_stages==2)] = "S_W"
  current_stages[which(current_stages==3)] = "I_W"
  current_stages[which(current_stages==4)] = "S_D"
  current_stages[which(current_stages==5)] = "I_D"
  current_stages = as.matrix(current_stages)
  
  return(current_stages)
}
###INITIAL_BEETLES
#
#set up the initial number of beetles in the right trees
initial_beetles = function(default_params, stages, IC_beetles){
  env = environment()
  list2env(default_params,env)
  
  #initialize the total population vector to 0 (categorized by tree)
  pop0ByTrees = rep(0,N*Nbs)
  
  Oi = rep(0,N)
  Oi[which(stages=="I_D")] = IC_beetles
  
  for (j in 1:N){
    pop0ByTrees[(j-1)*Nbs+2] = Oi[j]
  }
  return(as.matrix(pop0ByTrees))
}

###CLUSTER_INF_TREES
#
#select all trees in the circle of radius r and center center_world
cluster_inf_trees = function(elms, centre_world, r){
  ## Function to select some trees from the data (around the "center" within a radius r)
  diff_x = elms$X-centre_world$x
  diff_y = elms$Y-centre_world$y
  rad_xy = diff_x^2+diff_y^2
  idx_selected = which(rad_xy <= (r^2))
  return(idx_selected)
}

###PROBA_INFECTION
#
#Set the probability that, if a tree is susceptible and with beetles that are supposed to infect these trees, compute the proba that these trees are going to be infected next year (year 1)
#Remark: this situation is never used in the simulations of the paper so the output is always a vector full of zeros.
proba_infection = function(default_params, pop0ByTrees, stages){
  indexJi = seq(6,default_params$Nbs*default_params$N,by=default_params$Nbs)
  idx_presence_Ji = which(pop0ByTrees[indexJi]>0)
  
  vec.inf = rep(0,default_params$N)
  idx_H_or_Ws = which(stages=="H"|stages=="S_W")
  
  vecvec = intersect(idx_presence_Ji,idx_H_or_Ws)
  for (k in vecvec){
    vec.inf[k] = rbinom(1,pop0ByTrees[k],default.params$p_i)
    if(vec.inf[k]>0){
      vec.inf[k]=1
    }
  }
  return(vec.inf)
}

###CREATE_IC
#
#Function that creates the initial condition required with all info in an list of outputs
create_IC = function(sim_constants, 
                     elms,
                     IC_type,
                     IC_beetles,
                     IC_radius,
                     IC_number_dead_trees){
  
  centre_world = find_center(elms)
  # centre_world = special_center(elms)
  
  if (IC_type == "cluster"){ # just one cluster in the middle of the map
    sampling.dead.trees = cluster_inf_trees(elms,centre_world,IC_radius)
  } else if(IC_type == "2clusters"){ # 2 clusters
    center1 = special_center(elms) #center of the first cluster
    center2 = special_center2(elms) #center of the second cluster
    radius1 = IC_radius$r1
    radius2 = IC_radius$r2
    sampling.dead.trees1 = cluster_inf_trees(elms,center1,radius1)
    sampling.dead.trees2 = cluster_inf_trees(elms,center2,radius2)
    sampling.dead.trees = c(sampling.dead.trees1,sampling.dead.trees2)
  }else if(IC_type == "random"){ # dead trees are chosen randomly
    sampling.dead.trees = sample.int(sim_constants$default_params$N,size = IC_number_dead_trees)
  }
  
  IC_number_dead_trees = length(sampling.dead.trees) #this is the number of initially infected trees
  
  stages = initially_inf_trees(sim_constants$default_params,sampling.dead.trees) #IC for trees
  pop0ByTrees = initial_beetles(sim_constants$default_params,stages,IC_beetles) #IC for beetles
  infection0 = proba_infection(sim_constants$default_params,pop0ByTrees,stages) #infection for next year
  
  list.ic = list(stages=stages,pop0ByTrees=pop0ByTrees,IC_type=IC_type,IC_beetles=IC_beetles,IC_radius=IC_radius,IC_number_dead_trees=IC_number_dead_trees,infection0=infection0)
  return(list.ic)
}

### WEEK_TYPES
#
# Given a vector v of week numbers, return a vector of same size with the type of week under consideration
week_types = function(v, climate = NULL) {
  out = mat.or.vec(nr = length(v), nc = 1)
  weekType = list() #here is the regular year
  weekType[["Winter"]]        = c(1:21,45:53)
  weekType[["Emergence"]]     = 22
  weekType[["Breeding"]]      = 23:38
  weekType[["New_generation"]]     = 39:44
  for (wt in names(weekType)) {
    out[which(v %in% weekType[[wt]])] = wt
  }
  return(out)
}

### READ_PARAMETERS 
#
# read parameters from the csv file and replace those choosen in the main file (input)
# read the pre_processing files required for the given value of R_B (3 files for each R_B)
# read the tree database
# read the initial condition
# read the proba_root files created by the pre-processing of roots
# also adjust some parameters in function of the read files (nb of trees for instance)
read_parameters = function(sim_constants, input) {
  sim_constants$FILES = list()
  sim_constants$FILES[[1]] = sprintf("%s/parameters.csv", sim_constants$DIRS$DATA)
  parameters <- read.csv(sim_constants$FILES[[1]], header=TRUE)
  #### Sublist with parameters
  sim_constants$default_params = list()
  for (i in 1:dim(parameters)[1]) {
    sim_constants$default_params[[sprintf("%s",parameters[i,1])]] = 
      as.numeric(parameters[i,2])
  }
  # Override some values
  sim_constants$default_params$p_i  = input$p_i
  sim_constants$default_params$p_r  = input$p_r
  sim_constants$default_params$s_dt = input$s_dt
  sim_constants$default_params$R_B  = input$R_B
  
  ## Now we can set file names for the corresponding preprocessing since we have set the value of R_B
  sim_constants$FILES[[2]] = sprintf("%s/neighbours_%s_RB%d_%s.Rds", 
                                     sim_constants$DIRS$preproc_dists, 
                                     sim_constants$sim_core, 
                                     sim_constants$default_params$R_B,
                                     sim_constants$neighbourhood)
  sim_constants$FILES[[3]] = sprintf("%s/distance_neighbours_%s_RB%s_%s.Rds", 
                                     sim_constants$DIRS$preproc_dists, 
                                     sim_constants$sim_core, 
                                     sim_constants$default_params$R_B,
                                     sim_constants$neighbourhood)
  sim_constants$FILES[[4]] = sprintf("%s/neighbours_pos_%s_RB%s_%s.Rds", 
                                     sim_constants$DIRS$preproc_dists, 
                                     sim_constants$sim_core, 
                                     sim_constants$default_params$R_B,
                                     sim_constants$neighbourhood)
  # Also set file names for elms and roots files
  sim_constants$FILES[[5]] = sprintf("%s/elms_%s.Rds", 
                                     sim_constants$DIRS$nbhd_and_date, 
                                     sim_constants$neighbourhood)
  sim_constants$FILES[[7]] = sprintf("%s/proba_roots_%s.Rds",
                                     sim_constants$DIRS$nbhd_and_date,
                                     sim_constants$neighbourhood)
  
  # parameters <- read.csv(sim_constants$FILES[[1]], header=TRUE)
  neighbours_circle = readRDS(sim_constants$FILES[[2]])
  distance_neighbours = readRDS(sim_constants$FILES[[3]])
  neighbours_pos = readRDS(sim_constants$FILES[[4]])
  elms = readRDS(sim_constants$FILES[[5]])
  sim_constants$default_params$N = dim(elms)[1]
  sim_constants$default_params$proba_roots = readRDS(sim_constants$FILES[[7]])
  
  #### Sublist with preprocessed values
  sim_constants$default_params$preproc= list(neighbours_circle = neighbours_circle,
                                             distance_neighbours = distance_neighbours,
                                             neighbours_pos = neighbours_pos
  )
  #### Sublist with initial conditions
  # sim_constants$default_params$IC = readRDS(sim_constants$FILES[[6]])
  #### Sublist with selected elms
  sim_constants$default_params$elms = elms
  
  return(sim_constants)
}

### SET_OUTPUT_LOCATION
#
# set up the directory for the output with a name sim_date for the directory
# set_output_location = function(sim_constants) {
#   current_date_time = format(Sys.time(), "%Y_%m_%d_%H_%M")
#   abb_nbhd = abbreviate(sim_constants$Neighbourhood,minlength = 5)
#   # Create the folder for the output
#   output_dir = paste(sim_constants$DIRS$RESULTS,
#                      "/sim",
#                      current_date_time,
#                      "_",
#                      abb_nbhd,
#                      sep = "")
#   dir.create(output_dir)
#   sim_constants$DIRS$output_dir = output_dir
#   return(sim_constants)
# }

# SET_SIM_TIME
set_sim_time = function(sim_constants, startDate = start_date, endDate = end_date) {
  out = list()
  start_date = ISOweek2date(sprintf("%s-1",substr(date2ISOweek(startDate),1,8)))
  end_date =  ISOweek2date(sprintf("%s-1",substr(date2ISOweek(endDate),1,8)))  
  dates = readRDS(sprintf("%s/all_days_1958_to_2049.Rds",sim_constants$DIRS$DATA))
  weeksSims = seq(start_date, end_date, by = "week")
  idx = which(dates$dateFull %in% weeksSims)
  out = data.frame(idx = 1:length(weeksSims), dates[idx,])
  out$phase = week_types(out$weekShort)
  out$simple_year = as.numeric(out$year) - as.numeric(out$year[1])
  out$start_date = start_date
  out$end_date = end_date
  return(out)
}

# SET_SIM_ENVIRONMENT
# Returns a list with all the required parameters for a simulation, for each of the nb_sims simulations performed
set_sim_environment = function(sim_constants, input) {
  out = list()
  out$params = list()
  
  for (i in 1:sim_constants$nb_sims){
    out$params[[i]] = sim_constants$default_params
    varying_params = list()
    out$params[[i]]$varying_params = varying_params
    new_params = out$params[[i]]
    out$params[[i]]$matrices = set_demography_matrices(new_params)
  }
  return(out)
}

# SET_OTHER_CONSTANTS
#
# Gives constants that can be used and that will never be changed
# indexes and state names are defined here
set_other_constants = function(sim_constants) {
  indexOs = seq(from = 1,
                to = sim_constants$default_params$Nbs*sim_constants$default_params$N,
                by=sim_constants$default_params$Nbs)
  indexOi = seq(from = 2,
                to = sim_constants$default_params$Nbs*sim_constants$default_params$N,
                by=sim_constants$default_params$Nbs)
  indexMbs = seq(from = 3,
                 to = sim_constants$default_params$Nbs*sim_constants$default_params$N,
                 by=sim_constants$default_params$Nbs)
  indexMbi = seq(from = 4,
                 to = sim_constants$default_params$Nbs*sim_constants$default_params$N,
                 by=sim_constants$default_params$Nbs)
  indexJs = seq(from = 5,
                to = sim_constants$default_params$Nbs*sim_constants$default_params$N,
                by=sim_constants$default_params$Nbs)
  indexJi = seq(from = 6,
                to = sim_constants$default_params$Nbs*sim_constants$default_params$N,
                by=sim_constants$default_params$Nbs)
  indexMs = seq(from = 7,
                to = sim_constants$default_params$Nbs*sim_constants$default_params$N,
                by=sim_constants$default_params$Nbs)
  indexMi = seq(from = 8,
                to = sim_constants$default_params$Nbs*sim_constants$default_params$N,
                by=sim_constants$default_params$Nbs)
  indexAs = seq(from = 9,
                to = sim_constants$default_params$Nbs*sim_constants$default_params$N,
                by=sim_constants$default_params$Nbs)
  indexAi = seq(from = 10,
                to = sim_constants$default_params$Nbs*sim_constants$default_params$N,
                by=sim_constants$default_params$Nbs)
  statesNames <- c("H", "S_W", "I_W", "S_D", "I_D")
  
  out = list(statesNames = statesNames,
             indexOs = indexOs,
             indexOi = indexOi,
             indexMbs = indexMbs,
             indexJs = indexJs,
             indexMs = indexMs,
             indexAs = indexAs,
             indexMbi = indexMbi,
             indexJi = indexJi,
             indexMi = indexMi,
             indexAi = indexAi)
  return(out)
}

# SET_DEMOGRAPHY_MATRICES
#
# Define the demography matrices as a function of the tree status and the period
set_demography_matrices = function(default_params){
  env = environment()
  list2env(default_params,env)
  
  ##########################################################################
  ##matrices for winter event
  #######################################################################
  lH1 = c(0,0,0,0,0,0,0,0,0,0)
  lH2 = c(0,0,0,0,0,0,0,0,0,0)
  lH3 = c(0,0,0,0,0,0,0,0,0,0)
  lH4 = c(0,0,0,0,0,0,0,0,0,0)
  lH5 = c(0,0,0,0,S_J_winter*s_dt,0,0,0,0,0)
  lH6 = c(0,0,0,0,0,S_J_winter*s_dt,0,0,0,0)
  lH7 = c(0,0,0,0,0,0,0,0,0,0)
  lH8 = c(0,0,0,0,0,0,0,0,0,0)
  lH9 = c(0,0,0,0,0,0,0,0,0,0)
  lH10 = c(0,0,0,0,0,0,0,0,0,0)
  LH_winter = matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),
                     nrow = Nbs, ncol = Nbs, byrow = TRUE)
  LH_winter = as(LH_winter,"sparseMatrix")
  
  LSD_winter = mat.or.vec(Nbs,Nbs)
  LSD_winter = as(LSD_winter,"sparseMatrix")
  LID_winter = LSD_winter #beetles die if they overwinter in a dead or weak tree
  LIW_winter = LID_winter
  LSW_winter = LID_winter
  
  list_winter = list(LH_winter = LH_winter,
                     LSD_winter = LSD_winter,
                     LID_winter = LID_winter,
                     LSW_winter = LSW_winter,
                     LIW_winter = LIW_winter)
  
  #######################################################################
  ##matrices for Emergence event
  #######################################################################
  lH1 = c(0,0,0,0,0,0,0,0,0,0)
  lH2 = c(0,0,0,0,0,0,0,0,0,0)
  lH3 = c(0,0,0,0,0,0,0,0,0,0)
  lH4 = c(0,0,0,0,0,0,0,0,0,0)
  lH5 = c(0,0,0,0,s_dt*s_C,0,0,0,0,0)
  lH6 = c(0,0,0,0,0,s_dt*s_C,0,0,0,0)
  lH7 = c(0,0,0,0,s_dt*s_MC,0,0,0,0,0)
  lH8 = c(0,0,0,0,0,s_dt*s_MC,0,0,0,0)
  lH9 = c(0,0,0,0,0,0,0,0,0,0)
  lH10 = c(0,0,0,0,0,0,0,0,0,0)
  LH_emergence = matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),
                        nrow = Nbs, ncol = Nbs, byrow = TRUE)
  LH_emergence = as(LH_emergence,"sparseMatrix")
  
  lW1 = c(0,0,0,0,0,0,0,0,0,0)
  lW2 = c(0,0,0,0,0,0,0,0,0,0)
  lW3 = c(0,0,0,0,0,0,0,0,0,0)
  lW4 = c(0,0,0,0,0,0,0,0,0,0)
  lW5 = c(0,0,0,0,s_dt*s_C,0,0,0,0,0)
  lW6 = c(0,0,0,0,0,s_dt*s_C,0,0,0,0)
  lW7 = c(0,0,0,0,s_dt*s_MC,0,0,0,0,0)
  lW8 = c(0,0,0,0,0,s_dt*s_MC,0,0,0,0)
  lW9 = c(0,0,0,0,0,0,0,0,0,0)
  lW10 = c(0,0,0,0,0,0,0,0,0,0)
  LIW_emergence = matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),
                         nrow = Nbs, ncol = Nbs, byrow = TRUE)
  LIW_emergence = as(LIW_emergence, "sparseMatrix")
  LSW_emergence = LIW_emergence
  
  LSD_emergence = mat.or.vec(Nbs,Nbs)
  LSD_emergence = as(LSD_emergence,"sparseMatrix")
  LID_emergence = LSD_emergence #no beetle can be present in dead trees at this time of the year
  
  list_emergence = list(LH_emergence = LH_emergence,
                        LSW_emergence = LSW_emergence,
                        LIW_emergence = LIW_emergence,
                        LSD_emergence = LSD_emergence,
                        LID_emergence = LID_emergence)
  ###############################################################
  ##Matrices for breeding
  ###############################################################
  ####OLD CODE: les 2 lignes en commentaires etaient utilises avant de dire qu'on peut
  # lH1 = c(0,0,0,0,0,0,0,0,0,0)
  # lH2 = c(0,0,0,0,0,0,0,0,0,0)
  lH1 = c(s_dt*(s_J+s_JJp+s_FJ),0,0,0,0,0,0,0,f_JA,f_JA)
  lH2 = c(0,s_dt*(s_Jp+s_FJp),0,0,0,0,0,0,0,0)
  lH3 = c(0,0,0,0,0,0,0,0,0,0)
  lH4 = c(0,0,0,0,0,0,0,0,0,0)
  lH5 = c(0,0,0,0,0,0,0,0,0,0) # I remove juveniles
  lH6 = c(0,0,0,0,0,0,0,0,0,0)
  lH7 = c(0,0,0,0,s_dt,0,0,0,0,0) # I force juveniles to become ms or mi
  lH8 = c(0,0,0,0,0,s_dt,0,0,0,0)
  lH9 = c(0,0,0,0,0,0,0,0,0,0)
  lH10 = c(0,0,0,0,0,0,0,0,0,0)
  LH_breeding = matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),
                       nrow = Nbs, ncol = Nbs, byrow = TRUE)
  LH_breeding = as(LH_breeding,"sparseMatrix")
  
  lW1 = c(s_dt*(s_J+s_FJ),0,0,0,0,0,0,0,f_JA,f_JA)
  lW2 = c(s_dt*s_JJp,s_dt*(s_Jp+s_FJp),0,0,0,0,0,0,0,0)
  lW3 = c(0,0,0,0,0,0,0,0,0,0)# these four stages are not possible at this moment 
  lW4 = c(0,0,0,0,0,0,0,0,0,0)
  lW5 = c(0,0,0,0,0,0,0,0,0,0) 
  lW6 = c(0,0,0,0,0,0,0,0,0,0)
  lW7 = c(0,0,0,0,s_dt,0,0,0,0,0)
  lW8 = c(0,0,0,0,0,s_dt,0,0,0,0)
  lW9 = c(0,0,0,0,0,0,s_dt,0,0,0)
  lW10 = c(0,0,0,0,0,0,0,s_dt,0,0)
  LIW_breeding = matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),
                        nrow = Nbs, ncol = Nbs, byrow = TRUE)
  ## Check if this is good
  # LIW_breeding=LIW_breeding
  LIW_breeding = as(LIW_breeding,"sparseMatrix")
  
  lW1 = c(s_dt*(s_J+s_JJp+s_FJ),0,0,0,0,0,0,0,f_JA,f_JA)
  lW2 = c(0,s_dt*(s_Jp+s_FJp),0,0,0,0,0,0,0,0)
  lW3 = c(0,0,0,0,0,0,0,0,0,0)# these four stages are not possible at this moment
  lW4 = c(0,0,0,0,0,0,0,0,0,0)
  lW5 = c(0,0,0,0,0,0,0,0,0,0) 
  lW6 = c(0,0,0,0,0,0,0,0,0,0)
  lW7 = c(0,0,0,0,s_dt,0,0,0,0,0)
  lW8 = c(0,0,0,0,0,s_dt,0,0,0,0)
  lW9 = c(0,0,0,0,0,0,s_dt,0,0,0)
  lW10 = c(0,0,0,0,0,0,0,s_dt,0,0)
  LSW_breeding = matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),
                        nrow = Nbs, ncol = Nbs, byrow = TRUE)
  LSW_breeding = as(LSW_breeding,"sparseMatrix")
  #the same thing but s_JJp=0 (row2,col1)
  
  lD1 = c(0,0,0,0,0,0,0,0,f_JA,f_JA)
  lD2 = c(s_dt*(s_JJp+s_J+s_FJ),s_dt*(s_Jp+s_FJp),0,0,0,0,0,0,0,0)
  lD3 = c(0,0,0,0,0,0,0,0,0,0) #if the beetle Os is in an infected dead tree, then it becomes automatically an Oi (not an MBOs)
  lD4 = c(0,0,0,0,0,0,0,0,0,0) # not possible at this moment
  lD5 = c(0,0,0,0,0,0,0,0,0,0)
  lD6 = c(0,0,0,0,0,0,0,0,0,0)
  lD7 = c(0,0,0,0,0,0,0,0,0,0)
  lD8 = c(0,0,0,0,0,0,0,0,0,0)
  lD9 = c(0,0,0,0,0,0,s_dt,0,0,0)
  lD10 = c(0,0,0,0,0,0,0,s_dt,0,0)
  LID_breeding=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),
                      nrow=Nbs, ncol=Nbs, byrow = TRUE)
  LID_breeding=as(LID_breeding,"sparseMatrix")
  
  lD1 = c(s_dt*(s_J+s_FJ+s_JJp),0,0,0,0,0,0,0,f_JA,f_JA)
  lD2 = c(0,s_dt*(s_Jp+s_FJp),0,0,0,0,0,0,0,0)
  lD3 = c(0,0,0,0,0,0,0,0,0,0)
  lD4 = c(0,0,0,0,0,0,0,0,0,0)
  lD5 = c(0,0,0,0,0,0,0,0,0,0)
  lD6 = c(0,0,0,0,0,0,0,0,0,0)
  lD7 = c(0,0,0,0,0,0,0,0,0,0)
  lD8 = c(0,0,0,0,0,0,0,0,0,0)
  lD9 = c(0,0,0,0,0,0,s_dt,0,0,0)
  lD10 = c(0,0,0,0,0,0,0,s_dt,0,0)
  LSD_breeding=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),
                      nrow=Nbs, ncol=Nbs, byrow = TRUE)
  LSD_breeding=as(LSD_breeding,"sparseMatrix")
  list_breeding = list(LH_breeding = LH_breeding,
                       LSW_breeding = LSW_breeding,
                       LIW_breeding = LIW_breeding,
                       LSD_breeding = LSD_breeding,
                       LID_breeding = LID_breeding)
  
  ################################################################
  ##Matrices for new generation event
  ################################################################
  
  lH1 = c(s_dt*(s_J+s_JJp),0,0,0,0,0,0,0,0,0)
  lH2 = c(0,s_dt*s_Jp,0,0,0,0,0,0,0,0)
  lH3 = c(s_dt*s_FJ,0,0,0,0,0,0,0,0,0)
  lH4 = c(0,s_dt*s_FJp,0,0,0,0,0,0,0,0)
  lH5 = c(0,0,s_dt*s_CF,0,s_dt*s_C,0,0,0,0,0)
  lH6 = c(0,0,0,s_dt*s_CF,0,s_dt*s_C,0,0,0,0)
  lH7 = c(0,0,0,0,0,0,0,0,0,0) # I force callow adults to become ms or mi
  lH8 = c(0,0,0,0,0,0,0,0,0,0)
  lH9 = c(0,0,0,0,0,0,0,0,0,0)
  lH10 = c(0,0,0,0,0,0,0,0,0,0)
  LH_new_generation=matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),
                           nrow = Nbs, ncol = Nbs, byrow = TRUE)
  LH_new_generation=as(LH_new_generation, "sparseMatrix")
  
  lW1 = c(s_dt*s_J,0,0,0,0,0,0,0,0,0)
  lW2 = c(s_dt*s_JJp,s_dt*s_Jp,0,0,0,0,0,0,0,0)
  lW3 = c(s_dt*s_FJ,0,0,0,0,0,0,0,0,0)
  lW4 = c(0,s_dt*s_FJp,0,0,0,0,0,0,0,0)
  lW5 = c(0,0,s_dt*s_CF,0,s_dt*s_C,0,0,0,0,0)
  lW6 = c(0,0,0,s_dt*s_CF,0,s_dt*s_C,0,0,0,0)
  lW7 = c(0,0,0,0,0,0,0,0,0,0)
  lW8 = c(0,0,0,0,0,0,0,0,0,0)
  lW9 = c(0,0,0,0,0,0,0,0,0,0)
  lW10 = c(0,0,0,0,0,0,0,0,0,0)
  LIW_new_generation = matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),
                              nrow = Nbs,ncol = Nbs,byrow = TRUE)
  ## Check if this is good
  LIW_new_generation = as(LIW_new_generation, "sparseMatrix")
  
  lW1 = c(s_dt*(s_J+s_JJp),0,0,0,0,0,0,0,0,0)
  lW2 = c(0,s_dt*s_Jp,0,0,0,0,0,0,0,0)
  lW3 = c(s_dt*s_FJ,0,0,0,0,0,0,0,0,0)
  lW4 = c(0,s_dt*s_FJp,0,0,0,0,0,0,0,0)
  lW5 = c(0,0,s_dt*s_CF,0,s_dt*s_C,0,0,0,0,0)
  lW6 = c(0,0,0,s_dt*s_CF,0,s_dt*s_C,0,0,0,0)
  lW7 = c(0,0,0,0,0,0,0,0,0,0)
  lW8 = c(0,0,0,0,0,0,0,0,0,0)
  lW9 = c(0,0,0,0,0,0,0,0,0,0)
  lW10 = c(0,0,0,0,0,0,0,0,0,0)
  LSW_new_generation = matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),
                              nrow = Nbs, ncol = Nbs, byrow = TRUE)
  LSW_new_generation = as(LSW_new_generation,"sparseMatrix")
  #the same thing but s_JJp=0 (row2,col1)
  
  lD1 = c(0,0,0,0,0,0,0,0,0,0)
  lD2 = c(s_dt*(s_JJp+s_J),s_dt*s_Jp,0,0,0,0,0,0,0,0)
  lD3 = c(0,0,0,0,0,0,0,0,0,0)
  lD4 = c(s_dt*s_FJ,s_dt*s_FJp,0,0,0,0,0,0,0,0)
  lD5 = c(0,0,s_dt*s_CF,0,s_dt*s_C,0,0,0,0,0)
  lD6 = c(0,0,0,s_dt*s_CF,0,s_dt*s_C,0,0,0,0)
  lD7 = c(0,0,0,0,0,0,0,0,0,0)
  lD8 = c(0,0,0,0,0,0,0,0,0,0)
  lD9 = c(0,0,0,0,0,0,0,0,0,0)
  lD10 = c(0,0,0,0,0,0,0,0,0,0)
  LID_new_generation=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),
                            nrow=Nbs, ncol=Nbs, byrow = TRUE)
  LID_new_generation=as(LID_new_generation, "sparseMatrix")
  
  lD1 = c(s_dt*(s_J+s_JJp),0,0,0,0,0,0,0,0,0)
  lD2 = c(0,s_dt*s_Jp,0,0,0,0,0,0,0,0)
  lD3 = c(s_dt*s_FJ,0,0,0,0,0,0,0,0,0)
  lD4 = c(0,s_dt*s_FJp,0,0,0,0,0,0,0,0)
  lD5 = c(0,0,0,0,0,0,0,0,0,0)
  lD6 = c(0,0,0,0,0,0,0,0,0,0)
  lD7 = c(0,0,0,0,0,0,0,0,0,0)
  lD8 = c(0,0,0,0,0,0,0,0,0,0)
  lD9 = c(0,0,0,0,0,0,0,0,0,0)
  lD10 = c(0,0,0,0,0,0,0,0,0,0)
  LSD_new_generation=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),
                            nrow=Nbs, ncol=Nbs, byrow = TRUE)
  LSD_new_generation=as(LSD_new_generation,"sparseMatrix")
  
  list_new_generation = list(LH_new_generation = LH_new_generation,
                             LIW_new_generation = LIW_new_generation,
                             LSW_new_generation = LSW_new_generation,
                             LSD_new_generation = LSD_new_generation,
                             LID_new_generation = LID_new_generation)
  
  out = list(list_new_generation = list_new_generation,
             list_winter = list_winter,
             list_emergence = list_emergence,
             list_breeding = list_breeding)
  return(out)
}
