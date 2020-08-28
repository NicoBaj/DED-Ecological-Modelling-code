##pre_IC.R:
#functions that permit to create the three types of initial conditions (cluster, 2clusters and random)

###FIND_CENTER
#
#center of the map
find_center = function(Elms){
  min_X = min(Elms$X)
  max_X = max(Elms$X)
  min_Y = min(Elms$Y)
  max_Y = max(Elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  return(center)
}

###SPECIAL_CENTER
#
#set the center of the third quadrant (bottom left)
special_center = function(Elms){
  min_X = min(Elms$X)
  max_X = max(Elms$X)
  min_Y = min(Elms$Y)
  max_Y = max(Elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  third_quarter_center = list(x=(center$x+min_X)/2,y=(center$y+min_Y)/2)
  return(third_quarter_center)
}

###SPECIAL_CENTER2
#
#set the right center of the map (middle right)
special_center2 = function(Elms){
  min_X = min(Elms$X)
  max_X = max(Elms$X)
  min_Y = min(Elms$Y)
  max_Y = max(Elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  middle_right_center = list(x=(center$x+max_X)/2,y=(max_Y+min_Y)/2)
  return(middle_right_center)
}

###SPECIAL_CENTER3
#
#set the center of the first quadrant (top right)
special_center3 = function(Elms){
  min_X = min(Elms$X)
  max_X = max(Elms$X)
  min_Y = min(Elms$Y)
  max_Y = max(Elms$Y)
  center = list(x=(max_X+min_X)/2,y=(max_Y+min_Y)/2)
  first_quadrant_center = list(x=(center$x+max_X)/2,y=(center$y+max_Y)/2)
  return(first_quadrant_center)
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

###INITIALLY_INFECTED_TREES
#
#Set up the stages of trees at the initial time
initially_inf_trees = function(default_params,dead_sampling){
  current_stages = rep(1,default_params$N) #1 is for healthy
  initially_dead_trees_idx = dead_sampling 
  current_stages[initially_dead_trees_idx]=5 #5 is for dead and infected tree
  
  current_stages[which(current_stages==1)] = "H"
  current_stages[which(current_stages==2)] = "Ws"
  current_stages[which(current_stages==3)] = "Wi"
  current_stages[which(current_stages==4)] = "Ds"
  current_stages[which(current_stages==5)] = "Di"
  current_stages = as.matrix(current_stages)
  
  return(current_stages)
}
###INITIAL_BEETLES
#
#set up the initial number of beetles in the right trees
initial_beetles = function(default_params,stages,IC_beetles){
  env = environment()
  list2env(default_params,env)
  
  #initialize the total population vector to 0 (categorized by tree)
  pop0ByTrees = rep(0,N*Nbs)
  
  Oi = rep(0,N)
  Oi[which(stages=="Di")] = IC_beetles
  
  for (j in 1:N){
    pop0ByTrees[(j-1)*Nbs+2] = Oi[j]
  }
  return(as.matrix(pop0ByTrees))
}

###CLUSTER_INF_TREES
#
#select all trees in the circle of radius r and center center_world
cluster_inf_trees = function(Elms,centre_world,r){
  ## Function to select some trees from the data (around the "center" within a radius r)
  diff_x = Elms$X-centre_world$x
  diff_y = Elms$Y-centre_world$y
  rad_xy = diff_x^2+diff_y^2
  idx_selected = which(rad_xy <= (r^2))
  return(idx_selected)
}

###PROBA_INFECTION
#
#Set the probability that, if a tree is susceptible and with beetles that are supposed to infect these trees, compute the proba that these trees are going to be infected next year (year 1)
#Remark: this situation is never used in the simulations of the paper so the output is always a vector full of zeros.
proba_infection = function(default_params,pop0ByTrees,stages){
  indexJi = seq(6,default_params$Nbs*default_params$N,by=default_params$Nbs)
  idx_presence_Ji = which(pop0ByTrees[indexJi]>0)
  
  vec.inf = rep(0,default_params$N)
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

###CREATE_IC
#
#Function that creates the initial condition required with all info in an list of outputs
create_IC = function(sim_constants,Elms,IC_type,IC_beetles,IC_radius,IC_number_dead_trees){
  
  centre_world = find_center(Elms)
  # centre_world = special_center(Elms)
  
  if (IC_type == "cluster"){ # just one cluster in the middle of the map
    sampling.dead.trees = cluster_inf_trees(Elms,centre_world,IC_radius)
  } else if(IC_type == "2clusters"){ # 2 clusters
    center1 = special_center(Elms) #center of the first cluster
    center2 = special_center2(Elms) #center of the second cluster
    radius1 = IC_radius$r1
    radius2 = IC_radius$r2
    sampling.dead.trees1 = cluster_inf_trees(Elms,center1,radius1)
    sampling.dead.trees2 = cluster_inf_trees(Elms,center2,radius2)
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