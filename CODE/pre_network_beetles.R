##pre_network_beetles.R
# File that produces 

library(sqldf)
library(igraph)
library(Matrix)
library(R.utils)
library(parallel)

##############################################

### DISTANCE
#
# gives the distance between two trees
distance = function(xi,yi,xj,yj){
  d = sqrt((xi-xj)^2+(yi-yj)^2)
  return(d)
}

###NEIGHBOURS
#
#return the list of neighbours for each tree and the list of distances (for one value of R_B)
#note 1: that we first restrict the list through a square around the tree and then remove the trees too far, this allows to loop over a smaller number of trees
#note 2: the list of neighbours gives the real ID (from Tree.ID in the dataset)
neighbours = function(R_B,Elms){
  distance_neighbours = list() #necessary ?
  neighbours_square = list()
  neighbours_circle = list()
  neighbours.X = list()
  neighbours.Y = list()
  for (i in 1:dim(Elms)[1]){
    tree_ID = Elms$Tree.ID[i]
    x0 = Elms$X[i]
    y0 = Elms$Y[i]
    list.inside.X = which(abs(x0-Elms$X)<=R_B)
    list.inside.Y = which(abs(y0-Elms$Y)<=R_B)
    list.inside.square = intersect(list.inside.X,list.inside.Y)
    #list.inside.square is the set of trees that are inside the square of center (x0,y0) and side 2R_B
    
    neighbours_square[[i]] = setdiff(Elms$Tree.ID[list.inside.square],tree_ID) #gives the IDs of trees that are inside the square, except the tree itself
    list.inside.square = list.inside.square[! list.inside.square %in% i] #again, this is just to remove the value i from the vector
    
    neighbours.X[[i]] = Elms$X[list.inside.square]
    neighbours.Y[[i]] = Elms$Y[list.inside.square]
    
    distance_neighbours[[i]] = neighbours_square[[i]] #just to initialize at the same length
    ll = c() #vectors of trees in the square but not in the 
    
    if(length(neighbours_square[[i]])>0){
      for (j in 1:length(neighbours_square[[i]])){ 
        distance_neighbours[[i]][j] = distance(x0,y0,neighbours.X[[i]][j],neighbours.Y[[i]][j])
        if(!is.na(distance_neighbours[[i]][j])){
          if(distance_neighbours[[i]][j]>R_B){
            ll = c(ll,j)
          }
        }
      }
    }
    if(length(ll)>0){
      distance_neighbours[[i]]=distance_neighbours[[i]][-ll] #remove the neighbours too far (in the list of distances)
      neighbours_square[[i]] = neighbours_square[[i]][-ll] #remove the neighbours too far (in the list of neighbours)
    }
  }
  
  neighbours_circle = neighbours_square #After the removal, the neighbours_square is now a neighbours_circle
  
  return(list(neighbours_circ=neighbours_circle,distance_neighbours=distance_neighbours))
}

###NEIGHB_POS
#
#return the list of neighbours in the neighbourhood (not the entire city)
neighb_pos = function(Elms,neighbours_circle){#position of neighbours (return a list of vectors)
  lookup_ID = Elms$Tree.ID
  lookup_ID = cbind(lookup_ID,1:length(lookup_ID))#trees are positioned from 1 to N, N be the number of trees in the neighbourhood
  colnames(lookup_ID) = c("TreeID","idx")
  
  nb_pos = list()
  for (i in 1:length(neighbours_circle)) {#loop to give the position of each tree in the neighbourhood
    if (length(neighbours_circle[[i]])>0) {
      nb_pos[[i]] = neighbours_circle[[i]]*0 #initialize the vector to 0
      for (j in 1:length(neighbours_circle[[i]])) {
        pos = which(lookup_ID[,"TreeID"]==neighbours_circle[[i]][j])
        nb_pos[[i]][j] = lookup_ID[pos,"idx"]
      }
    }
    else {
      nb_pos[[i]] = mat.or.vec(1,0) #put the zero value, it is the same than put an empty list
    }
  }
  return(nb_pos)
}

###############################################################################

# Make code location aware using the library rprojroot
if(!require("rprojroot")){
  install.packages("rprojroot")
  require("rprojroot")
}
# Where is the current file? Needs to be sourced
this_file <- rprojroot::thisfile()
# Set code directory
TOP_DIR_CODE <- dirname(this_file)
# Set top directory
TOP_DIR <- dirname(TOP_DIR_CODE)

dir_data <- sprintf("%s/DATA",TOP_DIR)
dir_save <- sprintf("%s/Preprocessing/new_pre_processing",dir_data)

name_ngh = "NORTH_RIVER_HEIGHTS"

Elms = readRDS(sprintf("%s/Elms_Neighbourhood/Elms_%s.Rds",dir_data,name_ngh))

list.R_B = list()
seq_R_B = seq(20,100,by=20) #list of values of R_B
for (i in 1:length(seq_R_B)){
  list.R_B[[i]] = seq_R_B[i]
}
RUN_PARALLEL = FALSE
if(RUN_PARALLEL){
  no_cores <- detectCores()
  # Initiate cluster
  tictoc::tic()
  cl <- makeCluster(no_cores)
  clusterExport(cl,
                c("Elms"
                  ,"distance"
                  ,"neighbours"
                  ,"neighb_pos"
                ),
                envir = .GlobalEnv)
  # Run computation
  outputs = parLapply(cl = cl, X = list.R_B, fun =  function(x) neighbours(x,Elms))
  # Stop cluster
  stopCluster(cl)
  timeLoading=tictoc::toc()
}else{
  outputs = lapply(X = list.R_B, FUN =  function(x) neighbours(x,Elms))
}

for (i in 1:length(outputs)){
  output=outputs[[i]]
  R_B=seq_R_B[i]
  print(sprintf("R_B = %s",R_B))
  neighbours_circle = output[[1]]
  distance_neighbours = output[[2]]
  neighbours_pos = neighb_pos(Elms,output[[1]])
  Nb_trees = length(neighbours_circle)

  saveRDS(neighbours_circle,sprintf("%s/neighbours_%sTrees_RB%s_%s.Rds",dir_save,Nb_trees,R_B,name_ngh))
  saveRDS(distance_neighbours,sprintf("%s/distance_neighbours_%sTrees_RB%s_%s.Rds",dir_save,Nb_trees,R_B,name_ngh))
  saveRDS(neighbours_pos,sprintf("%s/neighbours_pos_%sTrees_RB%s_%s.Rds",dir_save,Nb_trees,R_B,name_ngh))
}
