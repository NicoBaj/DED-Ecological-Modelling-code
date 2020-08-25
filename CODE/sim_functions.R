########## Functions #########

### SYSTEM.OVER.TIME
#
# The main function that runs one simulation
system.over.time=function(sim_param,sim_constants){
  
  #Set up the initial conditions:
  matPopByTrees = mat.or.vec(sim_constants$default_params$N*sim_constants$default_params$Nbs,length(sim_constants$time$idx))
  colnames(matPopByTrees) = sim_constants$time$dateFull
  matPopByTrees[,1] = sim_param$IC$pop0ByTrees
  matPopByTrees<-as(matPopByTrees,"sparseMatrix")
  status_trees = mat.or.vec(sim_constants$default_params$N,length(sim_constants$time$idx))
  status_trees[,1] = sim_param$IC$stages
  
  new_infection = mat.or.vec(sim_constants$default_params$N,length(sim_constants$time$idx))
  new_infection[,1] = sim_param$IC$infection0
  
  event=sim_constants$time$phase[1]
  current_event = event
  #Demog is the big block diagonal demography matrix
  Demog = demography.matrices(sim_constants,sim_param$params,status_trees[,1],event)
  #End of IC
  
  root_or_vec = mat.or.vec(length(sim_constants$time$idx),2)
  
  #plot the intial set up
  if(PLOT_SIM){
    plot(sim_constants$default_params$elms$X,
         sim_constants$default_params$elms$Y,
         col="green",xlab = "X", ylab = "Y", main = "Trees")
    Di_Wi = which(status_trees[,1]=="Di"|status_trees[,1]=="Wi")
    points(sim_constants$default_params$elms$X[Di_Wi],sim_constants$default_params$elms$Y[Di_Wi],col="red")
  }
  

  #loop over time
  for (idx in 2:length(sim_constants$time$idx)) {
    event = sim_constants$time$phase[idx]
    # demography changes if event is different from before
    if (event != current_event) {
      Demog = demography.matrices(sim_constants,sim_param$params,status_trees[,idx-1],event)
      current_event = event
    }
    
    ##First movement
    matPopByTrees[,idx] = matPopByTrees[,idx-1]
    proba_infection=sim_constants$default_params$proba_infection
    mvt.Mbs.and.inf = mvt.beetles(sim_constants,sim_param$params,matPopByTrees[sim_constants$other$indexMbs,idx-1],status_trees[,(idx-1)],"H",VARYING_PREPROC,"no_carrier",proba_infection)
    mvt.Mbi.and.inf = mvt.beetles(sim_constants,sim_param$params,matPopByTrees[sim_constants$other$indexMbi,idx-1],status_trees[,(idx-1)],"H",VARYING_PREPROC,"carrier",proba_infection)
    mvt.Ms.and.inf = mvt.beetles(sim_constants,sim_param$params,matPopByTrees[sim_constants$other$indexMs,idx-1],status_trees[,(idx-1)],"D",VARYING_PREPROC,"no_carrier",proba_infection)
    mvt.Mi.and.inf = mvt.beetles(sim_constants,sim_param$params,matPopByTrees[sim_constants$other$indexMi,idx-1],status_trees[,(idx-1)],"D",VARYING_PREPROC,"carrier",proba_infection)
    
    matPopByTrees[sim_constants$other$indexMbs,idx] = mvt.Mbs.and.inf$Movement
    matPopByTrees[sim_constants$other$indexMbi,idx] = mvt.Mbi.and.inf$Movement
    matPopByTrees[sim_constants$other$indexMs,idx] = mvt.Ms.and.inf$Movement
    matPopByTrees[sim_constants$other$indexMi,idx] = mvt.Mi.and.inf$Movement
    
    #only the new feeders can infect trees
    new_infection[,idx]=mvt.Mbi.and.inf$Infection
    
    ## Then demography
    #Matrix Demog is defined first in the initial conditions, it only changes when scenario, trees status or year condition change
    matPopByTrees[,idx] = Demog%*%matPopByTrees[,idx]
    
    ##Finally update tree status
    if(event == "Winter" & sim_constants$time$phase[idx+1] == "Emerge" & idx != length(sim_constants$time$idx)){
      print(sprintf("Update tree states %d",idx))
      if(sim_constants$time$simple_year[idx]==1){#if this is the the first time we update
        vec_idx = 1:idx
      }else{#for the other years
        vec_idx = (idx-sim_constants$default_params$l):idx
      }
      idx_susceptible = which(status_trees[,idx-1]=="H"|status_trees[,idx-1]=="Ws")
      vec_new_inf = c()
      for (k in idx_susceptible){
        sum_proba_inf = sum(new_infection[k,vec_idx])
        if(sum_proba_inf>0){
          vec_new_inf = c(vec_new_inf,k)
        }
      }
      print(sprintf("nb of new infected trees by beetles = %s",length(vec_new_inf)))
      status_trees_after_beetles = mat.or.vec(sim_constants$default_params$N,1)
      status_trees_after_beetles[vec_new_inf] = "Wi"
      
      root_or_vec[idx,1] = length(vec_new_inf)
      
      #for all other trees, we use the markov chain as before (transition for ageing)
      other_trees = setdiff(1:sim_constants$default_params$N,vec_new_inf)
      
      for (k in other_trees){
        mcB <- new("markovchain", states = sim_constants$other$statesNames,
                   transitionMatrix = new.transition.matrices(sim_param$params))
        status_trees_after_beetles[k] = markovchainSequence(n = 1, markovchain = mcB, t0 = status_trees[k,(idx-1)])
      }
      
      ###############################
      # ROOT INFECTION
      ###############################
      
      if (sim_constants$roots){
        idx_Susceptible_trees = which(status_trees[,(idx-1)]=="H"|status_trees[,(idx-1)]=="Ws"|status_trees[,(idx-1)]=="Ds")
        status_trees_after_root = mat.or.vec(sim_constants$default_params$N,1)
        nb_inf_roots = 0
        
        for (j in idx_Susceptible_trees){#look up over susceptible trees
          ##j is the index of the tree in the "neighbourhood"
          # print(sprintf("Susceptible tree = %s",j))
          #first the rows that have j as a connected tree in proba_roots and find the associated neighbours
          
          #row_j is the indices of neighbours of j
          row_j = which(sim_param$params$proba_roots$idx_i==j)
          #pos_neighbours is the real position in the system
          pos_neighbours = sim_param$params$proba_roots$idx_j[row_j]
          
          #find the position of neighbours that are infected
          inf_neighbours = pos_neighbours[which(status_trees[pos_neighbours,idx-1]=="Di")]
          # print(sprintf("Number of infected ngh: %s",length(inf_neighbours)))
          # update_status=FALSE
          if(length(inf_neighbours)>0){
            #j_infected gives the rows for which a j is infected
            # points(sim_constants$default_params$elms$X[j],sim_constants$default_params$elms$Y[j],col="blue")
            #save the rows at which we have j associated to an infected neighbour
            
            row_inf = row_j[which(status_trees[sim_param$params$proba_roots$idx_j[row_j],idx-1]=="Di")]
            if(length(row_inf)>0){
              
              ##PB:
              vec_proba_inf = 2*sim_param$params$proba_roots$proba[row_inf]*sim_param$params$p_r
              tirage_poisson_binomiale = rpoibin(1,vec_proba_inf)
              if(tirage_poisson_binomiale>0){
                # update_status = TRUE
                if(status_trees[j,idx-1]=="H"|status_trees[j,idx-1]=="Ws"){
                  status_trees_after_root[j]="Wi"
                  nb_inf_roots = nb_inf_roots+1
                }else if(status_trees[j,idx-1]=="Ds"){
                  status_trees_after_root[j]="Di"
                  nb_inf_roots = nb_inf_roots+1
                }
              }
            }
          }
        }
        print(sprintf("nb of new infected trees by roots = %s",nb_inf_roots))
        root_or_vec[idx,2] = length(nb_inf_roots)
        
        status_trees[,idx] = merge_updates(sim_constants,status_trees_after_root,status_trees_after_beetles)
      
        
      }else{
        status_trees[,idx] = status_trees_after_beetles
        Di_Wi = which(status_trees[,idx]=="Di"|status_trees[,idx]=="Wi")
      }
      
      if (PLOT_SIM){
        plot(sim_constants$default_params$elms$X,sim_constants$default_params$elms$Y,
             col="green",xlab = "X", ylab = "Y", main = "Trees")
        Di_Wi = which(status_trees[,idx]=="Di"|status_trees[,idx]=="Wi")
        points(sim_constants$default_params$elms$X[Di_Wi],sim_constants$default_params$elms$Y[Di_Wi],col="red")
      }
      
      ##since we change the status, we change the demography matrices
      Demog=demography.matrices(sim_constants,sim_param$params,status_trees[,idx],event)
    }else{
      status_trees[,idx]=status_trees[,(idx-1)]
    }
  }
  sim_output = list(status_trees = status_trees)
  # sim_output = list(matPopByTrees=matPopByTrees,status_trees = status_trees)
  out = list(sim_output = sim_output,sim_param = sim_param,root_or_vec=root_or_vec)
  return(out)
}

convert_state_to_number = function(M){
  res = mat.or.vec(dim(M)[1],dim(M)[2])
  res[which(M=="H")] = 1
  res[which(M=="Ws")] = 2
  res[which(M=="Wi")] = 3
  res[which(M=="Ds")] = 4
  res[which(M=="Di")] = 5
  return(res)
}

merge_updates = function(sim_constants,root,beetles){
  #root and beetles are vectors of outcomes from both stochastic processes
  out = mat.or.vec(sim_constants$default_params$N,1)
  root = convert_state_to_number(as.matrix(root))
  beetles = convert_state_to_number(as.matrix(beetles))
  #indexes that are the same in both new updates
  idx_same_status = which(root==beetles)
  if(length(idx_same_status)>0){
    out[idx_same_status] = root[idx_same_status]
  }
  
  #indexes where root infection impacted more trees
  idx_diff_status_root = which(root>beetles)
  if(length(idx_diff_status_root)>0){
    out[idx_diff_status_root] = root[idx_diff_status_root]
  }
  
  #indexes where beetle infection impacted more trees
  idx_diff_status_beetles = which(beetles>root)
  if(length(idx_diff_status_beetles)>0){
    out[idx_diff_status_beetles] = beetles[idx_diff_status_beetles]
  }
  
  idx_root3 = which(root==3)
  idx_beetles4 = which(beetles==4)
  out[intersect(idx_root3,idx_beetles4)] = 5
  
  out[which(out==1)] = "H"
  out[which(out==2)] = "Ws"
  out[which(out==3)] = "Wi"
  out[which(out==4)] = "Ds"
  out[which(out==5)] = "Di"
  return(out)
}

mvt.beetles = function(sim_constants,params,vec.of.beetles,status_trees,type,VARYING_PREPROC,status_beetles,proba_infection){
  # params list of para: sim_param$params
  # vec.of.beetles vector of beetles is the vector of Mbs, Mbi, Ms or Mi
  # type MUST BE "H" or "D"
  
  ##old  
  # neighbours_circle = params$neighbours_circle
  # neighbours_pos = params$neighbours_pos
  # distance_neighbours = params$distance_neighbours
  
  ##New
  neighbours_circle = params$preproc$neighbours_circle
  neighbours_pos = params$preproc$neighbours_pos
  distance_neighbours = params$preproc$distance_neighbours
  
  result = rep(0,sim_constants$default_params$N)
  Infection = rep(0,sim_constants$default_params$N)
  vec.of.beetles = round(vec.of.beetles)
  list.random = list() 
  list.pos = list()
  list.dist = list()
  for (i in 1:sim_constants$default_params$N){
    if(vec.of.beetles[i]>0){
      if (type == "H"){## then beetles are looking for Healthy or Weak trees
        id_H = which(status_trees[neighbours_pos[[i]]]=="H")
        id_Ws = which(status_trees[neighbours_pos[[i]]]=="Ws")
        id_Wi = which(status_trees[neighbours_pos[[i]]]=="Wi")
        id = c(id_H,id_Ws,id_Wi)
        len = length(neighbours_pos[[i]][id])
      }
      else if (type == "D"){## then the beetles are looking for a dead or weak tree
        id_Ds = which(status_trees[neighbours_pos[[i]]]=="Ds")
        id_Di = which(status_trees[neighbours_pos[[i]]]=="Wi")
        id_Ws = which(status_trees[neighbours_pos[[i]]]=="Ws")
        id_Wi = which(status_trees[neighbours_pos[[i]]]=="Wi")
        id_H  = which(status_trees[neighbours_pos[[i]]]=="H")
        id = c(id_H,id_Ws,id_Wi,id_Ds,id_Di)
        len = length(neighbours_pos[[i]][id])
      }
      else{
        print(sprintf("no tree, type = %s, Tree number %s",type,i))
        len=0
      }
      
      if (len>0){
        pos_of_H_in_neighbourhood = neighbours_pos[[i]][id] # return the position of each neighbour
        distance_of_H_in_neighbourhood = distance_neighbours[[i]][id]
        
        #list.random[[i]] gives the destination of the moving beetles (this is random)
        list.random[[i]] = sample(1:len,size=vec.of.beetles[i],replace = TRUE) 
        #list.pos[[i]] gives the position of the destination trees
        list.pos[[i]] = pos_of_H_in_neighbourhood[list.random[[i]]]
        
        #we compute the proba to survive the distance for the destination trees
        # if(VARYING_PREPROC){
        #   list.dist[[i]] = proba.distance(params$maxD,distance_of_H_in_neighbourhood[list.random[[i]]])
        # }else{
        list.dist[[i]] = proba.distance(sim_constants$default_params$maxD,distance_of_H_in_neighbourhood[list.random[[i]]])
        # print("list.dist[[i]] for the ith tree")
        # print(sprintf("%s",list.dist[[i]]))
        # }
        #we then dispatch the beetles in the neighbour trees
        nb = list.pos[[i]] #POSITIONS
        dtfr = as.data.frame(table(nb))#create a data frame whose nb of rows = nb of neighbours
        dtfr$pb = as.data.frame(table(list.dist[[i]]))[,1] #we associate the proba to reach the associate tree
        dtfr$Survivals = 0 #
        if(status_beetles=="carrier" & type =="H"){#destination H, Ws or Wi
          dtfr$New_infection = 0
        }
        for (j in 1:dim(dtfr)[1]){
          dtfr$Survivals[j] = rbinom(1,dtfr$Freq[j],as.numeric(as.character(dtfr$pb[j])))
          if (status_beetles=="carrier" & type =="H"){
            if(dtfr$Survivals[j]>0){
              dtfr$New_infection[j] = rbinom(1,dtfr$Survival[j],proba_infection)
              if(is.na(dtfr$New_infection[j])){
                print(dtfr$Survival[j])
                print(proba_infection)
              }
              if(dtfr$New_infection[j]>0){
                dtfr$New_infection[j]=1
              }  
            }else{
              dtfr$New_infection[j]=0
            }
          }
          
        }
        result[as.numeric(as.vector(dtfr$nb))] = result[as.numeric(as.vector(dtfr$nb))] + dtfr$Survivals
        if(length(dtfr$New_infection)>0){
          Infection[as.numeric(as.vector(dtfr$nb))] = dtfr$New_infection
        }else{
          Infection[as.numeric(as.vector(dtfr$nb))]=0
        }
        
      }
    }
  }
  Infection[which(Infection>0)]=1
  out=list(Movement=result,Infection=Infection)
  return(out)
}

proba.distance = function(maxD,distance){
  res = exp(-distance/maxD)
  return(res)
}


### DEMOGRAPHY.MATRICES
#
# Set the big demography matrices when the period changes (Winter, Emerge, Breeding or Offspring) or when the status of trees change
# stagesi gives the new status of trees
# event gives the new scenario
# sim_param = sim_param$params
demography.matrices = function(sim_constants,sim_param,stagesi,event){#function to create the big demography matrix when status_trees or scenario change
  
  if (event == "Winter"){
    LH = sim_param$matrices$list_winter_normal$LH_win_normal
    LWs = sim_param$matrices$list_winter_normal$LWs_win_normal
    LWi = sim_param$matrices$list_winter_normal$LWi_win_normal
    LDs = sim_param$matrices$list_winter_normal$LDs_win_normal
    LDi = sim_param$matrices$list_winter_normal$LDi_win_normal
  }else if (event == "Emerge"){
    LH = sim_param$matrices$list_emergence1_normal$LH_em1_normal
    LWs = sim_param$matrices$list_emergence1_normal$LWs_em1_normal
    LWi = sim_param$matrices$list_emergence1_normal$LWi_em1_normal
    LDs = sim_param$matrices$list_emergence1_normal$LDs_em1_normal
    LDi = sim_param$matrices$list_emergence1_normal$LDi_em1_normal
  }else if (event == "Breeding"){
    LH = sim_param$matrices$list_breeding_normal$LH_bre_normal
    LWs = sim_param$matrices$list_breeding_normal$LWs_bre_normal
    LWi = sim_param$matrices$list_breeding_normal$LWi_bre_normal
    LDs = sim_param$matrices$list_breeding_normal$LDs_bre_normal
    LDi = sim_param$matrices$list_breeding_normal$LDi_bre_normal
  }else if (event =="Offspring"){
    LH = sim_param$matrices$list_offsprings_normal$LH_off_normal
    LWs = sim_param$matrices$list_offsprings_normal$LWs_off_normal
    LWi = sim_param$matrices$list_offsprings_normal$LWi_off_normal
    LDs = sim_param$matrices$list_offsprings_normal$LDs_off_normal
    LDi = sim_param$matrices$list_offsprings_normal$LDi_off_normal
  }
  
  L.list = list()
  for (i in (1:as.numeric(sim_constants$default_params$N))){
    s=stagesi[i]
    if(s=="H"){
      L.list[[i]] = LH
    }
    else if(s=="Ws"){
      L.list[[i]] = LWs
    }
    else if(s=="Wi"){
      L.list[[i]] = LWi
    }
    else if(s=="Ds"){
      L.list[[i]] = LDs
    }
    else if(s=="Di"){
      L.list[[i]] = LDi
    }
  }
  D=bdiag(L.list)
  # D = as.matrix(bdiag(L.list))
  D<-as(D,"sparseMatrix")
  return(D)
}

new.transition.matrices = function(params){
  env = environment()
  list2env(params,env)
  lt1 = c(Ph,0,0,0,0)
  lt2 = c(1-Ph,Pws,0,0,0)
  lt3 = c(0,0,0,0,0)
  lt4 = c(0,1-Pws,0,1,0)
  lt5 = c(0,0,1,0,1)
  Lt=matrix(data = c(lt1,lt2,lt3,lt4,lt5),nrow=5,ncol=5,byrow = TRUE)
  
  return(t(Lt))
}