###############################################
## Spread of the Dutch Elm Disease
## Code that creates initial conditions
## Code to select some trees
###############################################

##############################################
## Functions
##############################################

#center of the map, roughly
find_center = function(Elms){
  min_X = min(Elms$X)
  max_X = max(Elms$X)
  min_Y = min(Elms$Y)
  max_Y = max(Elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  return(center)
}

#center of the third quadrant
special_center = function(Elms){
  min_X = min(Elms$X)
  max_X = max(Elms$X)
  min_Y = min(Elms$Y)
  max_Y = max(Elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  first_quarter_center = list(x=(center$x+min_X)/2,y=(center$y+min_Y)/2)
  return(first_quarter_center)
  
}

special_center2 = function(Elms){
  min_X = min(Elms$X)
  max_X = max(Elms$X)
  min_Y = min(Elms$Y)
  max_Y = max(Elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  third_quarter_center = list(x=(center$x+max_X)/2,y=(max_Y+min_Y)/2)
  return(third_quarter_center)
}

special_center3 = function(Elms){
  min_X = min(Elms$X)
  max_X = max(Elms$X)
  min_Y = min(Elms$Y)
  max_Y = max(Elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  third_quarter_center = list(x=(center$x+max_X)/2,y=(center$y+max_Y)/2)
  return(third_quarter_center)
  
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

initially.inf.trees = function(default_params,dead_sampling){
  current_stages = rep(1,default_params$N) # trees are healthy
  initially_dead_trees_idx = dead_sampling 
  current_stages[initially_dead_trees_idx]=5 #Di
  
  current_stages[which(current_stages==1)] = "H"
  current_stages[which(current_stages==2)] = "Ws"
  current_stages[which(current_stages==3)] = "Wi"
  current_stages[which(current_stages==4)] = "Ds"
  current_stages[which(current_stages==5)] = "Di"
  current_stages = as.matrix(current_stages)
  
  return(current_stages)
}

initial.beetles = function(default_params,stages,IC_beetles){
  env = environment()
  list2env(default_params,env)
  
  pop0ByTrees = rep(0,N*Nbs) # initial population
  Os = rep(0,N)
  Oi = rep(0,N)
  Mbs = rep(0,N)
  Mbi = rep(0,N)
  Js = rep(0,N)
  Ji = rep(0,N)
  Ms = rep(0,N)
  Mi = rep(0,N)
  As = rep(0,N)
  Ai = rep(0,N)
  
  stages = convert_state_to_number(stages)
  idx_nonHealthy_trees = which(stages>1)
  Oi[which(stages==5)] = IC_beetles
  
  for(j in 1:N){
    pop0ByTrees[((j-1)*10+1):(j*10)] = c(Os[j],Oi[j],Mbs[j],Mbi[j],Js[j],Ji[j],Ms[j],Mi[j],As[j],Ai[j])
  }
  return(as.matrix(pop0ByTrees))
}

cluster.inf.trees = function(Elms,centre_world,r){
  ## Function to select some trees from the data (around the "center" within a radius r)
  
  SELECT_DISK = TRUE
  # r = 1100
  if (SELECT_DISK) {
    diff_x = Elms$X-centre_world$x
    diff_y = Elms$Y-centre_world$y
    rad_xy = diff_x^2+diff_y^2
    idx_selected = which(rad_xy <= (r^2))
  } else {
    diff_x = abs(Elms$X-centre_world$x)
    diff_y = abs(Elms$Y-centre_world$y)
    idx_selected = intersect(which(diff_x <= r),which(diff_y <= r))
  }
  
  return(idx_selected)
}

proba.infection = function(default_params,pop0ByTrees,stages){
  indexJi = seq(6,default_params$Nbs*default_params$N,by=default_params$Nbs)
  idx_presence_Ji = which(pop0ByTrees[indexJi]>0)
  vec.inf = rep(0,default_params$N)
  # stages = convert_state_to_number(stages)
  idx_H_or_Ws = which(stages=="H"|stages=="Ws")
  vecvec = intersect(idx_presence_Ji,idx_H_or_Ws)
  for (k in vecvec){
    vec.inf[k] = rbinom(1,pop0ByTrees[k],default.params$proba_infection)
    if(vec.inf[k]>0){
      vec.inf[k]=1
    }
  }
  return(vec.inf)
}

find_radius = function(Elms,center1,center2,IC_radius,nb_inf_trees){
  dead1 = cluster.inf.trees(Elms,center1,IC_radius)
  dead2 = cluster.inf.trees(Elms,center2,IC_radius)
  l = length(dead1)+length(dead2)
  if(l!=nb_inf_trees){
    cpt = 1
    while (l!=nb_inf_trees & cpt<100) {
      if(l>nb_inf_trees){
        IC_radius = IC_radius-0.5
        dead1 = cluster.inf.trees(Elms,center1,IC_radius)
        dead2 = cluster.inf.trees(Elms,center2,IC_radius)
        l = length(dead1)+length(dead2)
      }else{
        IC_radius = IC_radius+0.5
        dead1 = cluster.inf.trees(Elms,center1,IC_radius)
        dead2 = cluster.inf.trees(Elms,center2,IC_radius)
        l = length(dead1)+length(dead2)
      }
      cpt = cpt+1
    }
    out = IC_radius
  }else{
    out = IC_radius
  }
  return(out)
}

create_IC = function(sim_constants,Elms,IC_type,IC_beetles,IC_radius,IC_number_dead_trees){
  
  centre_world = find_center(Elms)
  # centre_world = special_center(Elms)
  
  if (IC_type == "cluster"){ # just one cluster in the middle of the map
    sampling.dead.trees = cluster.inf.trees(Elms,centre_world,IC_radius)
  } else if(IC_type == "2clusters"){ # 2 clusters
    center1 = special_center(Elms) #center of the first cluster
    center2 = special_center2(Elms) #center of the second cluster
    radius1 = IC_radius$r1
    radius2 = IC_radius$r2
    sampling.dead.trees1 = cluster.inf.trees(Elms,center1,radius1)
    sampling.dead.trees2 = cluster.inf.trees(Elms,center2,radius2)
    sampling.dead.trees = c(sampling.dead.trees1,sampling.dead.trees2)
  }else if(IC_type == "random"){ # dead trees are chosen randomly
    sampling.dead.trees = sample.int(sim_constants$default_params$N,size = IC_number_dead_trees)
    
    # if(is.integer(IC_number_dead_trees)){
    #   sampling.dead.trees = sample.int(sim_constants$default_params$N,size = IC_number_dead_trees)
    # }else{#if IC_number_dead_trees == "ic_radius": we want to have the same number of dead trees than a cluster of size IC_radius
    #   # cluster virtually created from the IC_radius 
    #   virtual_cluster = cluster.inf.trees(Elms,centre_world,IC_radius)
    #   # then, the length of virtual_cluster gives the number of dead trees we want
    #   sampling.dead.trees = sample.int(sim_constants$default_params$N,size = length(virtual_cluster))
    # }
  }
  
  stages = initially.inf.trees(sim_constants$default_params,sampling.dead.trees)
  pop0ByTrees = initial.beetles(sim_constants$default_params,stages,IC_beetles) #beetles initially present
  infection0 = proba.infection(sim_constants$default_params,pop0ByTrees,stages)
  
  list.ic = list(stages=stages,pop0ByTrees=pop0ByTrees,IC_type=IC_type,IC_beetles=IC_beetles,IC_radius=IC_radius,IC_number_dead_trees=IC_number_dead_trees,infection0=infection0)
  return(list.ic)
}

## keep te following commented part, this is when we want to create specific ic:
# Elms = readRDS("/storage/var/groups/mathbio/DED_DATA_OUTPUT/DATA/Elms_Neighbourhood/Elms_MIXED_PULBERRY_CRESCENT_PARK.RData")

# ic = create_IC(sim_constants,Elms,IC_type,IC_beetles,IC_radius,IC_number_dead_trees)
# saveRDS(ic,file = sprintf("/storage/var/groups/mathbio/DED_DATA_OUTPUT/PARAMS/IC_%s_radius_%s.RData",sim_constants$Neighbourhood,IC_radius))
# if(IC_type=="2clusters"){
#   center1 = special_center(Elms)
#   center2 = special_center3(Elms)
#   new_radius = find_radius(Elms,center1,center2,IC_radius,IC_number_dead_trees)
# }
