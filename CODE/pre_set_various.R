### SET_DIRECTORIES for outputs and data
set_directories = function(DIR_TOP) {
  out = list()
  out$TOP = DIR_TOP
  out$DATA = sprintf("%s/DATA",DIR_TOP)
  out$RESULT = sprintf("%s/RESULTS",DIR_TOP)
  return(out)
}


### WEEK_TYPES
#
# Given a vector v of week numbers, return a vector of same size with the type of week under consideration
week_types = function(v,climate = NULL) {
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
read_parameters = function(sim_constants,input) {
  sim_constants$FILES = list()
  sim_constants$FILES[[1]] = sprintf("%s/parameters.csv", sim_constants$DIRS$DATA)
  parameters <- read.csv(sim_constants$FILES[[1]], header=TRUE)
  #### Sublist with parameters
  sim_constants$default_params = list()
  for (i in 1:dim(parameters)[1]) {
    sim_constants$default_params[[sprintf("%s",parameters[i,1])]] = as.numeric(parameters[i,2])
  }
  sim_constants$default_params$p_i             = input$p_i
  sim_constants$default_params$p_r             = input$p_r
  sim_constants$default_params$s_dt             = input$s_dt
  sim_constants$default_params$R_B             = input$R_B
  
  ## Now we can load the good preprocessing since we have the right R_B
  sim_constants$FILES[[2]] = sprintf("%s/Preprocessing/neighbours_%s_maxD%s.RData", sim_constants$DIRS$DATA, sim_constants$sim_core, sim_constants$default_params$R_B)
  sim_constants$FILES[[3]] = sprintf("%s/Preprocessing/distance_neighbours_%s_maxD%s.RData", sim_constants$DIRS$DATA, sim_constants$sim_core, sim_constants$default_params$R_B)
  sim_constants$FILES[[4]] = sprintf("%s/Preprocessing/neighbours_pos_%s_maxD%s.RData", sim_constants$DIRS$DATA,sim_constants$sim_core, sim_constants$default_params$R_B)
  sim_constants$FILES[[5]] = sprintf("%s/Elms_Neighbourhood/Elms_%s.RData", sim_constants$DIRS$DATA, sim_constants$Neighbourhood)
  
  #the following file is just a default_file, then it changes when we set up the IC
  # sim_constants$FILES[[6]] = sprintf("%s/IC_PULBERRY_radius_60.RData", sim_constants$DIRS$DATA)
  
  # sim_constants$FILES[[7]] = sprintf("%s/Proba_roots/Proba_roots_%s_pr%s.RData",sim_constants$DIRS$DATA,sim_constants$Neighbourhood,sim_constants$default_params$p_r*100)
  
  sim_constants$FILES[[7]] = sprintf("%s/Proba_roots/Proba_roots_%s.Rds",sim_constants$DIRS$DATA,sim_constants$Neighbourhood)
  
  # parameters <- read.csv(sim_constants$FILES[[1]], header=TRUE)
  neighbours_circle = readRDS(sim_constants$FILES[[2]])
  distance_neighbours = readRDS(sim_constants$FILES[[3]])
  neighbours_pos = readRDS(sim_constants$FILES[[4]])
  Elms = readRDS(sim_constants$FILES[[5]])
  sim_constants$default_params$N = dim(Elms)[1]
  sim_constants$default_params$proba_roots = readRDS(sim_constants$FILES[[7]])
  
  #### Sublist with preprocessed values
  sim_constants$default_params$preproc= list(neighbours_circle = neighbours_circle,
                                             distance_neighbours = distance_neighbours,
                                             neighbours_pos = neighbours_pos
  )
  
  #### Sublist with initial conditions
  # sim_constants$default_params$IC = readRDS(sim_constants$FILES[[6]])
  
  #### Sublist with selected elms
  sim_constants$default_params$elms = Elms
  
  return(sim_constants)
}

### SET_OUTPUT_LOCATION
#
# set up the directory for the output with a name sim_date for the repertory
set_output_location = function(sim_constants) {
  current_date_time = format(Sys.time(), "%Y_%m_%d_%H_%M")
  abb_ngh = abbreviate(sim_constants$Neighbourhood,minlength = 5)
  
  #we create the folder for the outputs
  output_dir = paste(sim_constants$DIRS$RESULT,
                     "/sim",
                     current_date_time,
                     "_",
                     abb_ngh,
                     sep = "")
  dir.create(output_dir)
  sim_constants$output_dir = output_dir
  return(sim_constants)
}

# SET_SIM_TIME
set_sim_time = function(sim_constants, startDate = start_date, endDate = end_date) {
  out = list()
  start_date = ISOweek2date(sprintf("%s-1",substr(date2ISOweek(startDate),1,8)))
  end_date =  ISOweek2date(sprintf("%s-1",substr(date2ISOweek(endDate),1,8)))  
  dates = readRDS(sprintf("%s/allDays_1958_to_2049.RData",sim_constants$DIRS$DATA))
  weeksSims = seq(start_date, end_date, by = "week")
  idx = which(dates$dateFull %in% weeksSims)
  out = data.frame(idx = 1:length(weeksSims), dates[idx,])
  out$phase=week_types(out$weekShort)
  out$simple_year=as.numeric(out$year) - as.numeric(out$year[1])
  out$start_date=start_date
  out$end_date=end_date
  return(out)
}

# SET_SIM_ENVIRONMENT
# Returns a list with all the required parameters for a simulation, for each of the nb_sims simulations performed
set_sim_environment = function(sim_constants,input) {
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
  
  indexOs = seq(1,sim_constants$default_params$Nbs*sim_constants$default_params$N,by=sim_constants$default_params$Nbs)
  indexOi = seq(2,sim_constants$default_params$Nbs*sim_constants$default_params$N,by=sim_constants$default_params$Nbs)
  indexMbs = seq(3,sim_constants$default_params$Nbs*sim_constants$default_params$N,by=sim_constants$default_params$Nbs)
  indexMbi = seq(4,sim_constants$default_params$Nbs*sim_constants$default_params$N,by=sim_constants$default_params$Nbs)
  indexJs = seq(5,sim_constants$default_params$Nbs*sim_constants$default_params$N,by=sim_constants$default_params$Nbs)
  indexJi = seq(6,sim_constants$default_params$Nbs*sim_constants$default_params$N,by=sim_constants$default_params$Nbs)
  indexMs = seq(7,sim_constants$default_params$Nbs*sim_constants$default_params$N,by=sim_constants$default_params$Nbs)
  indexMi = seq(8,sim_constants$default_params$Nbs*sim_constants$default_params$N,by=sim_constants$default_params$Nbs)
  indexAs = seq(9,sim_constants$default_params$Nbs*sim_constants$default_params$N,by=sim_constants$default_params$Nbs)
  indexAi = seq(10,sim_constants$default_params$Nbs*sim_constants$default_params$N,by=sim_constants$default_params$Nbs)
  statesNames <- c("H", "S_W", "Wi", "Ds", "Di")
  # scenarii = c("bad","good","normal")
  out = list(statesNames=statesNames,indexOs=indexOs,indexOi=indexOi,indexMbs=indexMbs,indexJs=indexJs,indexMs=indexMs,indexAs=indexAs,indexMbi=indexMbi,indexJi=indexJi,indexMi=indexMi,indexAi=indexAi)
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
  LH_winter=matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LH_winter=as(LH_winter,"sparseMatrix")
  
  LSD_winter = mat.or.vec(Nbs,Nbs)
  LSD_winter = as(LSD_winter,"sparseMatrix")
  LID_winter = LSD_winter #beetles die if they overwinter in a dead or weak tree
  LIW_winter = LID_winter
  LSW_winter = LID_winter
  
  list_winter = list(LH_winter=LH_winter,LSD_winter=LSD_winter,LID_winter=LID_winter,LSW_winter=LSW_winter,LIW_winter=LIW_winter)
  
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
  LH_emergence=matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LH_emergence=as(LH_emergence,"sparseMatrix")
  
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
  LIW_emergence=matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LIW_emergence=as(LIW_emergence,"sparseMatrix")
  LSW_emergence=LIW_emergence
  
  LSD_emergence=mat.or.vec(Nbs,Nbs)
  LSD_emergence=as(LSD_emergence,"sparseMatrix")
  LID_emergence = LSD_emergence #no beetle can be present in dead trees at this time of the year
  
  list_emergence = list(LH_emergence=LH_emergence,LSW_emergence=LSW_emergence,LIW_emergence=LIW_emergence,LSD_emergence=LSD_emergence,LID_emergence=LID_emergence)
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
  LH_breeding=matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LH_breeding=as(LH_breeding,"sparseMatrix")
  
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
  LIW_breeding=matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  ## Check if this is good
  # LIW_breeding=LIW_breeding
  LIW_breeding=as(LIW_breeding,"sparseMatrix")
  
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
  LSW_breeding=matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LSW_breeding=as(LSW_breeding,"sparseMatrix")
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
  LID_breeding=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
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
  LSD_breeding=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LSD_breeding=as(LSD_breeding,"sparseMatrix")
  
  list_breeding = list(LH_breeding=LH_breeding,LSW_breeding=LSW_breeding,LIW_breeding=LIW_breeding,LSD_breeding=LSD_breeding,LID_breeding=LID_breeding)
  
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
  LH_new_generation=matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LH_new_generation=as(LH_new_generation,"sparseMatrix")
  
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
  LIW_new_generation=matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  ## Check if this is good
  LIW_new_generation=as(LIW_new_generation,"sparseMatrix")
  
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
  LSW_new_generation=matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LSW_new_generation=as(LSW_new_generation,"sparseMatrix")
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
  LID_new_generation=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LID_new_generation=as(LID_new_generation,"sparseMatrix")
  
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
  LSD_new_generation=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LSD_new_generation=as(LSD_new_generation,"sparseMatrix")
  
  list_new_generation = list(LH_new_generation=LH_new_generation,LIW_new_generation=LIW_new_generation,LSW_new_generation=LSW_new_generation,LSD_new_generation=LSD_new_generation,LID_new_generation=LID_new_generation)
  
  out = list(list_new_generation=list_new_generation,list_winter=list_winter,list_emergence=list_emergence,list_breeding=list_breeding
             )
  return(out)
}
