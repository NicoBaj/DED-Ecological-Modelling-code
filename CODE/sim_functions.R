########## Functions #########

### SYSTEM.OVER.TIME
#
# Runs one simulation
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
  Demog = demography_matrix(sim_constants,sim_param$params,status_trees[,1],event)
  #End of set up
  
  root_or_beetle = mat.or.vec(length(sim_constants$time$idx),2)
  
  #plot the intial set up
  if(sim_constants$GATES$PLOT_SIM){
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
      Demog = demography_matrix(sim_constants,sim_param$params,status_trees[,idx-1],event)
      current_event = event
    }
    
    ##First movement
    matPopByTrees[,idx] = matPopByTrees[,idx-1]
    # proba_infection = sim_constants$default_params$proba_infection
    mvt.Mbs.and.inf = mvt_beetles(sim_constants,sim_param$params,matPopByTrees[sim_constants$other$indexMbs,idx-1],status_trees[,(idx-1)],"H",VARYING_PREPROC,"no_carrier",sim_constants$default_params$proba_infection)
    mvt.Mbi.and.inf = mvt_beetles(sim_constants,sim_param$params,matPopByTrees[sim_constants$other$indexMbi,idx-1],status_trees[,(idx-1)],"H",VARYING_PREPROC,"carrier",sim_constants$default_params$proba_infection)
    mvt.Ms.and.inf = mvt_beetles(sim_constants,sim_param$params,matPopByTrees[sim_constants$other$indexMs,idx-1],status_trees[,(idx-1)],"all",VARYING_PREPROC,"no_carrier",sim_constants$default_params$proba_infection)
    mvt.Mi.and.inf = mvt_beetles(sim_constants,sim_param$params,matPopByTrees[sim_constants$other$indexMi,idx-1],status_trees[,(idx-1)],"all",VARYING_PREPROC,"carrier",sim_constants$default_params$proba_infection)
    
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
    if(event == "Winter" & sim_constants$time$phase[idx+1] == "Emergence" & idx != length(sim_constants$time$idx)){
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
      
      root_or_beetle[idx,1] = length(vec_new_inf)
      
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
        root_or_beetle[idx,2] = length(nb_inf_roots)
        
        status_trees[,idx] = merge_updates(sim_constants,status_trees_after_root,status_trees_after_beetles)
        
        
      }else{
        status_trees[,idx] = status_trees_after_beetles
        Di_Wi = which(status_trees[,idx]=="Di"|status_trees[,idx]=="Wi")
      }
      
      if (sim_constants$GATES$PLOT_SIM){
        plot(sim_constants$default_params$elms$X,sim_constants$default_params$elms$Y,
             col="green",xlab = "X", ylab = "Y", main = "Trees")
        Di_Wi = which(status_trees[,idx]=="Di"|status_trees[,idx]=="Wi")
        points(sim_constants$default_params$elms$X[Di_Wi],sim_constants$default_params$elms$Y[Di_Wi],col="red")
      }
      
      ##since we change the status, we change the demography matrices
      Demog=demography_matrix(sim_constants,sim_param$params,status_trees[,idx],event)
    }else{
      status_trees[,idx]=status_trees[,(idx-1)]
    }
  }
  sim_output = list(status_trees = status_trees)
  # sim_output = list(matPopByTrees=matPopByTrees,status_trees = status_trees)
  out = list(sim_output = sim_output,sim_param = sim_param,root_or_beetle=root_or_beetle)
  return(out)
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

###MVT_BEETLES
#
#function that moves beetles from trees to other. If beetles are carrying the spores, then the probability to infect the destination trees are computed and saved
mvt_beetles = function(sim_constants,params,vec.of.beetles,status_trees,type,VARYING_PREPROC,status_beetles,proba_infection){
  # params list of para: sim_param$params
  # vec.of.beetles vector of beetles is the vector of Mbs, Mbi, Ms or Mi
  # type MUST BE "H" or "all": this is the type of trees towards which the beetles go, "H" meaning not dead trees (H, SW and IW)
  
  #pre_processing from beetles are required, the following is just to reduce the names when call:
  neighbours_circle = params$preproc$neighbours_circle
  neighbours_pos = params$preproc$neighbours_pos
  distance_neighbours = params$preproc$distance_neighbours
  
  #The entry i of result gives the number of beetles that moved in tree i from all other possible trees 
  result = rep(0,sim_constants$default_params$N)
  #The entry i of Infection is either 0 or 1, if 1 then the tree i has been infected by a beetle, 0 otherwise
  Infection = rep(0,sim_constants$default_params$N)
  vec.of.beetles = round(vec.of.beetles)#
  list.random = list() 
  list.pos = list()
  list.dist = list()
  
  for (i in 1:sim_constants$default_params$N){
    if(vec.of.beetles[i]>0){#if beetles need to move :
      if (type == "H"){## then beetles are looking for Healthy or Weak trees
        id_H = which(status_trees[neighbours_pos[[i]]]=="H")
        id_Ws = which(status_trees[neighbours_pos[[i]]]=="Ws")
        id_Wi = which(status_trees[neighbours_pos[[i]]]=="Wi")
        id = c(id_H,id_Ws,id_Wi)
        len = length(neighbours_pos[[i]][id])
      }
      else if (type == "all"){## then the beetles are looking for any type of tree
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
      #len is the number of trees that beetles in tree i can reach
      #id is their index
      
      if (len>0){#if trees are reachable
        pos_of_H_in_neighbourhood = neighbours_pos[[i]][id] #position of each neighbour
        distance_of_H_in_neighbourhood = distance_neighbours[[i]][id] #distance of each neighbour
        
        list.random[[i]] = sample(1:len,size=vec.of.beetles[i],replace = TRUE)#destination for the beetles choosen randomly
        list.pos[[i]] = pos_of_H_in_neighbourhood[list.random[[i]]]#position of the destination trees
        
        list.dist[[i]] = proba.distance(sim_constants$default_params$maxD,distance_of_H_in_neighbourhood[list.random[[i]]]) #probability to survive the distance for the destination trees

        #then, dispatch the beetles in the neighbour trees
        nb = list.pos[[i]] 
        dtfr = as.data.frame(table(nb))#create a data frame whose nb of rows = nb of neighbours
        dtfr$pb = as.data.frame(table(list.dist[[i]]))[,1] #associate the proba to reach the associate tree
        dtfr$Survivals = 0 #initialize the beetles that survive to the travel to 0
        if(status_beetles=="carrier" & type =="H"){#destination H, Ws or Wi AND beetles are carrying the spores
          dtfr$New_infection = 0 #initialize the vector of trees that are going to be infected to 0
        }
        
        for (j in 1:dim(dtfr)[1]){
          dtfr$Survivals[j] = rbinom(1,dtfr$Freq[j],as.numeric(as.character(dtfr$pb[j])))
          if (status_beetles=="carrier" & type =="H"){#only carrier beetles can infect susceptible trees
            if(dtfr$Survivals[j]>0){#only beetles who survive can infect the destination tree ...
              dtfr$New_infection[j] = rbinom(1,dtfr$Survival[j],proba_infection)
              if(dtfr$New_infection[j]>0){#one successful infection suffices to infect the tree
                dtfr$New_infection[j]=1
              }  
            }else{
              dtfr$New_infection[j]=0
            }
          }
        }
        
        #each destination tree receives the beetles that survived the travel:
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
  #If a tree i is successfully infected by a beetle, then Infection[i]=1, else Infection[i]=0. This will be used to update the tree status at the right time
  out=list(Movement=result,Infection=Infection)
  return(out)
}

proba.distance = function(maxD,distance){
  res = exp(-distance/maxD)
  return(res)
}

### DEMOGRAPHY_MATRIX
#
# Set the big demography matrix when the period changes (Winter, Emerge, Breeding or Offspring) or when the status of trees change
# stagesi gives the new status of trees
# event gives the new scenario
# sim_param = sim_param$params
demography_matrix = function(sim_constants,sim_param,stagesi,event){#function to create the big demography matrix when status_trees or scenario change
  
  if (event == "Winter"){
    LH = sim_param$matrices$list_winter$LH_winter
    LSW = sim_param$matrices$list_winter$LSW_winter
    LIW = sim_param$matrices$list_winter$LID_winter
    LSD = sim_param$matrices$list_winter$LSD_winter
    LID = sim_param$matrices$list_winter$LID_winter
  }else if (event == "Emergence"){
    LH = sim_param$matrices$list_emergence$LH_emergence
    LSW = sim_param$matrices$list_emergence$LSW_emergence
    LIW = sim_param$matrices$list_emergence$LIW_emergence
    LSD = sim_param$matrices$list_emergence$LSD_emergence
    LID = sim_param$matrices$list_emergence$LID_emergence
  }else if (event == "Breeding"){
    LH = sim_param$matrices$list_breeding$LH_breeding
    LSW = sim_param$matrices$list_breeding$LSW_breeding
    LIW = sim_param$matrices$list_breeding$LIW_breeding
    LSD = sim_param$matrices$list_breeding$LSD_breeding
    LID = sim_param$matrices$list_breeding$LID_breeding
  }else if (event =="New_generation"){
    LH = sim_param$matrices$list_new_generation$LH_new_generation
    LSW = sim_param$matrices$list_new_generation$LSW_new_generation
    LIW = sim_param$matrices$list_new_generation$LIW_new_generation
    LSD = sim_param$matrices$list_new_generation$LSD_new_generation
    LID = sim_param$matrices$list_new_generation$LID_new_generation
  }
  
  L.list = list()
  for(i in 1:sim_constants$default_params$N){
    s=stagesi[i]
    if(s=="H"){
      L.list[[i]] = LH
    }
    else if(s=="Ws"){
      L.list[[i]] = LSW
    }
    else if(s=="Wi"){
      L.list[[i]] = LIW
    }
    else if(s=="Ds"){
      L.list[[i]] = LSD
    }
    else if(s=="Di"){
      L.list[[i]] = LID
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