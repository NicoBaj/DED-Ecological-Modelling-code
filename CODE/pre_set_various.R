### SET_DIRECTORIES for outputs and data
set_directories = function(DIR_TOP) {
  out = list()
  out$TOP = DIR_TOP
  out$DATA = sprintf("%s/DATA",DIR_TOP)
  out$PARAM = sprintf("%s/PARAMS",DIR_TOP)
  out$RESULT = sprintf("%s/RESULTS",DIR_TOP)
  return(out)
}


### WEEK_TYPES
#
# For the moment, we assume the schedule of weeks is fixed in the year. In the future, the schedule might vary.
# Given a vector v of week numbers, return a vector of same size with the type of week under consideration
week_types = function(v,climate = NULL) {
  out = mat.or.vec(nr = length(v), nc = 1)
  weekType = list() #here is the regular year
  weekType[["Winter"]]        = c(1:21,45:53)
  weekType[["Emerge"]]        = 22
  weekType[["Breeding"]]      = 23:38
  weekType[["Offspring"]]     = 39:44
  for (wt in names(weekType)) {
    out[which(v %in% weekType[[wt]])] = wt
  }
  return(out)
}

### READ_PARAMETERS
read_parameters = function(sim_constants,input) {
  sim_constants$FILES = list()
  sim_constants$FILES[[1]] = sprintf("%s/parameters.csv", sim_constants$DIRS$PARAM)
  parameters <- read.csv(sim_constants$FILES[[1]], header=TRUE)
  #### Sublist with parameters
  sim_constants$default_params = list()
  for (i in 1:dim(parameters)[1]) {
    sim_constants$default_params[[sprintf("%s",parameters[i,1])]] = as.numeric(parameters[i,2])
  }
  sim_constants$default_params$proba_infection = input$pi
  sim_constants$default_params$p_r             = input$pr
  sim_constants$default_params$Sdt             = input$sdt
  sim_constants$default_params$maxD            = input$maxD
  ## Now we can load the good preprocessing since we have the right maxD
  sim_constants$FILES[[2]] = sprintf("%s/Preprocessing/neighbours_%s_maxD%s.RData", sim_constants$DIRS$DATA, sim_constants$sim_core, sim_constants$default_params$maxD)
  sim_constants$FILES[[3]] = sprintf("%s/Preprocessing/distance_neighbours_%s_maxD%s.RData", sim_constants$DIRS$DATA, sim_constants$sim_core, sim_constants$default_params$maxD)
  sim_constants$FILES[[4]] = sprintf("%s/Preprocessing/neighbours_pos_%s_maxD%s.RData", sim_constants$DIRS$DATA,sim_constants$sim_core, sim_constants$default_params$maxD)
  sim_constants$FILES[[5]] = sprintf("%s/Elms_Neighbourhood/Elms_%s.RData", sim_constants$DIRS$DATA, sim_constants$Neighbourhood)
  #the following file is just a default_file, then it changes when we set up the IC
  sim_constants$FILES[[6]] = sprintf("%s/IC_PULBERRY_radius_60.RData", sim_constants$DIRS$PARAM)
  sim_constants$FILES[[7]] = sprintf("%s/Proba_roots/Proba_roots_%s_pr%s.RData",sim_constants$DIRS$DATA,sim_constants$Neighbourhood,sim_constants$default_params$p_r*100)
  
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
  sim_constants$default_params$IC = readRDS(sim_constants$FILES[[6]])
  #### Sublist with selected elms
  sim_constants$default_params$elms = Elms
  
  return(sim_constants)
}

### SET_OUTPUT_LOCATION
set_output_location = function(sim_constants) {
  current_date_time = format(Sys.time(), "%Y_%m_%d_%H_%M")
  # if (sim_constants$roots){
  #we create the folder for the outputs
  #   output_dir = sprintf(sprintf("%s/sim_%s_maxD%03g_pb_inf%03g_pr%03g_%s",sim_constants$DIRS$RESULT,sim_constants$initial_set_up$IC_type,sim_constants$default_params$maxD,sim_constants$default_params$proba_infection*100,sim_constants$default_params$p_r*100,sim_constants$sim_core))
  # }else{
  #   output_dir = sprintf(sprintf("%s/sim_%s_maxD%03g_pb_inf%03g_pr0_%s",sim_constants$DIRS$RESULT,sim_constants$initial_set_up$IC_type,sim_constants$default_params$maxD,sim_constants$default_params$proba_infection*100,sim_constants$sim_core))
  # }
  
  #we create the folder for the outputs
  output_dir = paste(sim_constants$DIRS$RESULT,
                     "/sim",
                     current_date_time,
                     "_",
                     sim_constants$sim_core,
                     sep = "")
  dir.create(output_dir)
  dir.create(paste(output_dir,"/imgs",sep=""))
  dir.create(paste(output_dir,"/proportions",sep=""))
  sim_constants[["filename"]] = paste(output_dir,"/sim.RData",sep="")
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
  #######################################################################################
  ## Provide the range and the name of the parameters -> the order matters !!!!!!!!!!!!!!
  # save the hypercube in the output file
  
  out$params[[1]] = sim_constants$default_params
  
      # Changing small parameters
      varying_params = list()
      out$params[[1]]$varying_params = varying_params
      new_params = out$params[[1]]
      out$params[[1]]$matrices = set_demography_matrices(new_params)
  return(out)
}

proba_infection = function(sim_constants,fileNb){
  nb.of.files = 10
  vec.pb = seq(0.01,0.1,length.out = nb.of.files)
  sim_constants$default_params$proba_infection = vec.pb[as.numeric(fileNb)]
  return(vec.pb[as.numeric(fileNb)])
}

set_maxD = function(sim_constants,fileNb){
  nb.of.files = 10
  # nb.of.files = 3
  vec.maxD = seq(20,380,by = 40)
  # vec.maxD = c(20,180,340)
  sim_constants$default_params$maxD = vec.maxD[as.numeric(fileNb)]
  return(vec.maxD[as.numeric(fileNb)])
}

set_maxD_IC_random = function(sim_constants,fileNb){
  nb.of.files = 3
  vec.maxD = c(300,300,300)
  maxD = vec.maxD[as.numeric(fileNb)]
  return(maxD)
}

# SET_OTHER_CONSTANTS
# Gives constants that can be used and that will never be changed
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
  statesNames <- c("H", "Ws", "Wi", "Ds", "Di")
  # scenarii = c("bad","good","normal")
  out = list(statesNames=statesNames,indexOs=indexOs,indexOi=indexOi,indexMbs=indexMbs,indexJs=indexJs,indexMs=indexMs,indexAs=indexAs,indexMbi=indexMbi,indexJi=indexJi,indexMi=indexMi,indexAi=indexAi)
  return(out)
}

set_demography_matrices = function(default_params){
  env = environment()
  list2env(default_params,env)
  ##########################################################################
  ##matrices for winter event, mean year
  #######################################################################
  lH1 = c(0,0,0,0,0,0,0,0,0,0)
  lH2 = c(0,0,0,0,0,0,0,0,0,0)
  lH3 = c(0,0,0,0,0,0,0,0,0,0)
  lH4 = c(0,0,0,0,0,0,0,0,0,0)
  lH5 = c(0,0,0,0,Sj_winter*Sdt,0,0,0,0,0)
  lH6 = c(0,0,0,0,0,Sj_winter*Sdt,0,0,0,0)
  lH7 = c(0,0,0,0,0,0,0,0,0,0)
  lH8 = c(0,0,0,0,0,0,0,0,0,0)
  lH9 = c(0,0,0,0,0,0,0,0,0,0)
  lH10 = c(0,0,0,0,0,0,0,0,0,0)
  LH_win_normal=matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LH_win_normal=as(LH_win_normal,"sparseMatrix")
  
  LDs_win_normal = mat.or.vec(Nbs,Nbs)
  LDs_win_normal = as(LDs_win_normal,"sparseMatrix")
  LDi_win_normal = LDs_win_normal #beetles die if they overwinter in a dead or weak tree
  LWi_win_normal = LDi_win_normal
  LWs_win_normal = LDi_win_normal
  
  list_winter_normal = list(LH_win_normal=LH_win_normal,LDs_win_normal=LDs_win_normal,LDi_win_normal=LDi_win_normal,LWs_win_normal=LWs_win_normal,LWi_win_normal=LWi_win_normal)
  
  #######################################################################
  ##matrices for Emergence1 event, mean year
  #######################################################################
  lH1 = c(0,0,0,0,0,0,0,0,0,0)
  lH2 = c(0,0,0,0,0,0,0,0,0,0)
  lH3 = c(0,0,0,0,0,0,0,0,0,0)
  lH4 = c(0,0,0,0,0,0,0,0,0,0)
  lH5 = c(0,0,0,0,Sdt*Sjs,0,0,0,0,0)
  lH6 = c(0,0,0,0,0,Sdt*Sji,0,0,0,0)
  lH7 = c(0,0,0,0,Sdt*SMJs,0,Sdt*Sms,0,0,0)
  lH8 = c(0,0,0,0,0,Sdt*SMJi,0,Sdt*Smi,0,0)
  lH9 = c(0,0,0,0,0,0,0,0,0,0)
  lH10 = c(0,0,0,0,0,0,0,0,0,0)
  LH_em1_normal=matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LH_em1_normal=as(LH_em1_normal,"sparseMatrix")
  
  lW1 = c(0,0,0,0,0,0,0,0,0,0)
  lW2 = c(0,0,0,0,0,0,0,0,0,0)
  lW3 = c(0,0,0,0,0,0,0,0,0,0)
  lW4 = c(0,0,0,0,0,0,0,0,0,0)
  lW5 = c(0,0,0,0,Sdt*Sjs,0,0,0,0,0)
  lW6 = c(0,0,0,0,0,Sdt*Sji,0,0,0,0)
  lW7 = c(0,0,0,0,Sdt*SMJs,0,Sdt*Sms,0,0,0)
  lW8 = c(0,0,0,0,0,Sdt*SMJi,0,Sdt*Smi,0,0)
  lW9 = c(0,0,0,0,0,0,0,0,0,0)
  lW10 = c(0,0,0,0,0,0,0,0,0,0)
  LWi_em1_normal=matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LWi_em1_normal=as(LWi_em1_normal,"sparseMatrix")
  LWs_em1_normal=LWi_em1_normal
  
  LDs_em1_normal=mat.or.vec(Nbs,Nbs)
  LDs_em1_normal=as(LDs_em1_normal,"sparseMatrix")
  LDi_em1_normal = LDs_em1_normal #no beetle can be present in dead trees at this time of the year
  
  list_emergence1_normal = list(LH_em1_normal=LH_em1_normal,LWs_em1_normal=LWs_em1_normal,LWi_em1_normal=LWi_em1_normal,LDs_em1_normal=LDs_em1_normal,LDi_em1_normal=LDi_em1_normal)
  ###############################################################
  ##Matrices for breeding
  ###############################################################
  ####OLD CODE: les 2 lignes en commentaires etaient utilises avant de dire qu'on peut
  # lH1 = c(0,0,0,0,0,0,0,0,0,0)
  # lH2 = c(0,0,0,0,0,0,0,0,0,0)
  lH1 = c(Sdt*(Sos+Sosoi+SmbOs),0,0,0,0,0,0,0,Fecs,Fecs)
  lH2 = c(0,Sdt*(Soi+SmbOi),0,0,0,0,0,0,0,0)
  lH3 = c(0,0,0,0,0,0,0,0,0,0)
  lH4 = c(0,0,0,0,0,0,0,0,0,0)
  lH5 = c(0,0,0,0,0,0,0,0,0,0) # I remove juveniles
  lH6 = c(0,0,0,0,0,0,0,0,0,0)
  lH7 = c(0,0,0,0,Sdt,0,Sdt*Sms,0,0,0) # I force juveniles to become ms or mi
  lH8 = c(0,0,0,0,0,Sdt,0,Sdt*Smi,0,0)
  lH9 = c(0,0,0,0,0,0,0,0,0,0)
  lH10 = c(0,0,0,0,0,0,0,0,0,0)
  LH_bre_normal=matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LH_bre_normal=as(LH_bre_normal,"sparseMatrix")
  
  lW1 = c(Sdt*(Sos+SmbOs),0,0,0,0,0,0,0,Fecs,Fecs)
  lW2 = c(Sdt*Sosoi,Sdt*(Soi+SmbOi),0,0,0,0,0,0,0,0)
  lW3 = c(0,0,0,0,0,0,0,0,0,0)# these four stages are not possible at this moment 
  lW4 = c(0,0,0,0,0,0,0,0,0,0)
  lW5 = c(0,0,0,0,0,0,0,0,0,0) 
  lW6 = c(0,0,0,0,0,0,0,0,0,0)
  lW7 = c(0,0,0,0,Sdt,0,Sdt*Sms,0,0,0)
  lW8 = c(0,0,0,0,0,Sdt,0,Sdt*Smi,0,0)
  lW9 = c(0,0,0,0,0,0,Sdt*SAMs,0,Sdt*Sas,0)
  lW10 = c(0,0,0,0,0,0,0,Sdt*SAMi,0,Sdt*Sai)
  LWi_bre_normal=matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  ## Check if this is good
  # LWi_bre_normal=LWi_bre_normal
  LWi_bre_normal=as(LWi_bre_normal,"sparseMatrix")
  
  lW1 = c(Sdt*(Sos+Sosoi+SmbOs),0,0,0,0,0,0,0,Fecs,Fecs)
  lW2 = c(0,Sdt*(Soi+SmbOi),0,0,0,0,0,0,0,0)
  lW3 = c(0,0,0,0,0,0,0,0,0,0)# these four stages are not possible at this moment
  lW4 = c(0,0,0,0,0,0,0,0,0,0)
  lW5 = c(0,0,0,0,0,0,0,0,0,0) 
  lW6 = c(0,0,0,0,0,0,0,0,0,0)
  lW7 = c(0,0,0,0,Sdt,0,Sdt*Sms,0,0,0)
  lW8 = c(0,0,0,0,0,Sdt,0,Sdt*Smi,0,0)
  lW9 = c(0,0,0,0,0,0,Sdt*SAMs,0,Sdt*Sas,0)
  lW10 = c(0,0,0,0,0,0,0,Sdt*SAMi,0,Sdt*Sai)
  LWs_bre_normal=matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LWs_bre_normal=as(LWs_bre_normal,"sparseMatrix")
  #the same thing but SOsOi=0 (row2,col1)
  
  lD1 = c(0,0,0,0,0,0,0,0,Fecs,Fecs)
  lD2 = c(Sdt*(Sosoi+Sos+SmbOs),Sdt*(Soi+SmbOi),0,0,0,0,0,0,0,0)
  lD3 = c(0,0,0,0,0,0,0,0,0,0) #if the beetle Os is in an infected dead tree, then it becomes automatically an Oi (not an MBOs)
  lD4 = c(0,0,0,0,0,0,0,0,0,0) # not possible at this moment
  lD5 = c(0,0,0,0,0,0,0,0,0,0)
  lD6 = c(0,0,0,0,0,0,0,0,0,0)
  lD7 = c(0,0,0,0,0,0,0,0,0,0)
  lD8 = c(0,0,0,0,0,0,0,0,0,0)
  lD9 = c(0,0,0,0,0,0,Sdt*SAMs,0,Sdt*Sas,0)
  lD10 = c(0,0,0,0,0,0,0,Sdt*SAMi,0,Sdt*Sai)
  LDi_bre_normal=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LDi_bre_normal=as(LDi_bre_normal,"sparseMatrix")
  
  lD1 = c(Sdt*(Sos+SmbOs+Sosoi),0,0,0,0,0,0,0,Fecs,Fecs)
  lD2 = c(0,Sdt*(Soi+SmbOi),0,0,0,0,0,0,0,0)
  lD3 = c(0,0,0,0,0,0,0,0,0,0)
  lD4 = c(0,0,0,0,0,0,0,0,0,0)
  lD5 = c(0,0,0,0,0,0,0,0,0,0)
  lD6 = c(0,0,0,0,0,0,0,0,0,0)
  lD7 = c(0,0,0,0,0,0,0,0,0,0)
  lD8 = c(0,0,0,0,0,0,0,0,0,0)
  lD9 = c(0,0,0,0,0,0,Sdt*SAMs,0,Sdt*Sas,0)
  lD10 = c(0,0,0,0,0,0,0,Sdt*SAMi,0,Sdt*Sai)
  LDs_bre_normal=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LDs_bre_normal=as(LDs_bre_normal,"sparseMatrix")
  
  list_breeding_normal = list(LH_bre_normal=LH_bre_normal,LWs_bre_normal=LWs_bre_normal,LWi_bre_normal=LWi_bre_normal,LDs_bre_normal=LDs_bre_normal,LDi_bre_normal=LDi_bre_normal)
  
  ################################################################
  ##Matrices for Offsprings development
  ################################################################
  
  #######OLD CODE
  # lH1 = c(0,0,0,0,0,0,0,0,0,0)
  # lH2 = c(0,0,0,0,0,0,0,0,0,0)
  # lH3 = c(0,0,0,0,0,0,0,0,0,0)
  # lH4 = c(0,0,0,0,0,0,0,0,0,0)
  # lH5 = c(0,0,Sdt*SJMbs,0,Sdt*Sjs,0,0,0,0,0)
  # lH6 = c(0,0,0,Sdt*SJMbi,0,Sdt*Sji,0,0,0,0)
  # lH7 = c(0,0,0,0,0,0,0,0,0,0) # I force callow adults to become ms or mi
  # lH8 = c(0,0,0,0,0,0,0,0,0,0)
  # lH9 = c(0,0,0,0,0,0,0,0,0,0)
  # lH10 = c(0,0,0,0,0,0,0,0,0,0)
  
  lH1 = c(Sdt*(Sos+Sosoi),0,0,0,0,0,0,0,0,0)
  lH2 = c(0,Sdt*Soi,0,0,0,0,0,0,0,0)
  lH3 = c(Sdt*SmbOs,0,Sdt*Smbs,0,0,0,0,0,0,0)
  lH4 = c(0,Sdt*SmbOi,0,Sdt*Smbi,0,0,0,0,0,0)
  lH5 = c(0,0,Sdt*SJMbs,0,Sdt*Sjs,0,0,0,0,0)
  lH6 = c(0,0,0,Sdt*SJMbi,0,Sdt*Sji,0,0,0,0)
  lH7 = c(0,0,0,0,0,0,0,0,0,0) # I force callow adults to become ms or mi
  lH8 = c(0,0,0,0,0,0,0,0,0,0)
  lH9 = c(0,0,0,0,0,0,0,0,0,0)
  lH10 = c(0,0,0,0,0,0,0,0,0,0)
  LH_off_normal=matrix(data = c(lH1,lH2,lH3,lH4,lH5,lH6,lH7,lH8,lH9,lH10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LH_off_normal=as(LH_off_normal,"sparseMatrix")
  
  lW1 = c(Sdt*Sos,0,0,0,0,0,0,0,0,0)
  lW2 = c(Sdt*Sosoi,Sdt*Soi,0,0,0,0,0,0,0,0)
  lW3 = c(Sdt*SmbOs,0,Sdt*Smbs,0,0,0,0,0,0,0)
  lW4 = c(0,Sdt*SmbOi,0,Sdt*Smbi,0,0,0,0,0,0)
  lW5 = c(0,0,Sdt*SJMbs,0,Sdt*Sjs,0,0,0,0,0)
  lW6 = c(0,0,0,Sdt*SJMbi,0,Sdt*Sji,0,0,0,0)
  lW7 = c(0,0,0,0,0,0,0,0,0,0)
  lW8 = c(0,0,0,0,0,0,0,0,0,0)
  lW9 = c(0,0,0,0,0,0,0,0,0,0)
  lW10 = c(0,0,0,0,0,0,0,0,0,0)
  LWi_off_normal=matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  ## Check if this is good
  LWi_off_normal=as(LWi_off_normal,"sparseMatrix")
  
  lW1 = c(Sdt*(Sos+Sosoi),0,0,0,0,0,0,0,0,0)
  lW2 = c(0,Sdt*Soi,0,0,0,0,0,0,0,0)
  lW3 = c(Sdt*SmbOs,0,Sdt*Smbs,0,0,0,0,0,0,0)
  lW4 = c(0,Sdt*SmbOi,0,Sdt*Smbi,0,0,0,0,0,0)
  lW5 = c(0,0,Sdt*SJMbs,0,Sdt*Sjs,0,0,0,0,0)
  lW6 = c(0,0,0,Sdt*SJMbi,0,Sdt*Sji,0,0,0,0)
  lW7 = c(0,0,0,0,0,0,0,0,0,0)
  lW8 = c(0,0,0,0,0,0,0,0,0,0)
  lW9 = c(0,0,0,0,0,0,0,0,0,0)
  lW10 = c(0,0,0,0,0,0,0,0,0,0)
  LWs_off_normal=matrix(data = c(lW1,lW2,lW3,lW4,lW5,lW6,lW7,lW8,lW9,lW10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LWs_off_normal=as(LWs_off_normal,"sparseMatrix")
  #the same thing but SOsOi=0 (row2,col1)
  
  lD1 = c(0,0,0,0,0,0,0,0,0,0)
  lD2 = c(Sdt*(Sosoi+Sos),Sdt*Soi,0,0,0,0,0,0,0,0)
  lD3 = c(0,0,Sdt*Smbs,0,0,0,0,0,0,0)
  lD4 = c(Sdt*SmbOs,Sdt*SmbOi,0,Sdt*Smbi,0,0,0,0,0,0)
  lD5 = c(0,0,Sdt*SJMbs,0,Sdt*Sjs,0,0,0,0,0)
  lD6 = c(0,0,0,Sdt*SJMbi,0,Sdt*Sji,0,0,0,0)
  lD7 = c(0,0,0,0,0,0,0,0,0,0)
  lD8 = c(0,0,0,0,0,0,0,0,0,0)
  lD9 = c(0,0,0,0,0,0,0,0,0,0)
  lD10 = c(0,0,0,0,0,0,0,0,0,0)
  LDi_off_normal=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LDi_off_normal=as(LDi_off_normal,"sparseMatrix")
  
  lD1 = c(Sdt*(Sos+Sosoi),0,0,0,0,0,0,0,0,0)
  lD2 = c(0,Sdt*Soi,0,0,0,0,0,0,0,0)
  lD3 = c(Sdt*SmbOs,0,Sdt*Smbs,0,0,0,0,0,0,0)
  lD4 = c(0,Sdt*SmbOi,0,Sdt*Smbi,0,0,0,0,0,0)
  lD5 = c(0,0,0,0,0,0,0,0,0,0)
  lD6 = c(0,0,0,0,0,0,0,0,0,0)
  lD7 = c(0,0,0,0,0,0,0,0,0,0)
  lD8 = c(0,0,0,0,0,0,0,0,0,0)
  lD9 = c(0,0,0,0,0,0,0,0,0,0)
  lD10 = c(0,0,0,0,0,0,0,0,0,0)
  LDs_off_normal=matrix(data = c(lD1,lD2,lD3,lD4,lD5,lD6,lD7,lD8,lD9,lD10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  LDs_off_normal=as(LDs_off_normal,"sparseMatrix")
  
  list_offsprings_normal = list(LH_off_normal=LH_off_normal,LWi_off_normal=LWi_off_normal,LWs_off_normal=LWs_off_normal,LDs_off_normal=LDs_off_normal,LDi_off_normal=LDi_off_normal)
  
  lC1 = c(0,0,0,0,0,0,0,0,0,0)
  lC2 = c(0,0,0,0,0,0,0,0,0,0)
  lC3 = c(0,0,0,0,0,0,0,0,0,0)
  lC4 = c(0,0,0,0,0,0,0,0,0,0)
  lC5 = c(0,0,0,0,0,0,0,0,0,0)
  lC6 = c(0,0,0,0,0,0,0,0,0,0)
  lC7 = c(0,0,0,0,0,0,0,0,0,0)
  lC8 = c(0,0,0,0,0,0,0,0,0,0)
  lC9 = c(0,0,0,0,0,0,0,0,0,0)
  lC10 = c(0,0,0,0,0,0,0,0,0,0)
  LC=matrix(data = c(lC1,lC2,lC3,lC4,lC5,lC6,lC7,lC8,lC9,lC10),nrow=Nbs,ncol=Nbs,byrow = TRUE)
  # LCT=t(LC)
  LC=as(LC,"sparseMatrix")
  
  out = list(list_offsprings_normal=list_offsprings_normal,list_winter_normal=list_winter_normal,list_emergence1_normal=list_emergence1_normal,list_breeding_normal=list_breeding_normal,LC=LC)
  return(out)
}
