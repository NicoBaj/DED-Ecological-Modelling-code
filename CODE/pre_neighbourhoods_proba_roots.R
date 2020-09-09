##pre_load_trees_and_separate_into_ngh.R
#load the right dataset
#split the dataset for each neighbourhood

####################################################################################
### SELECTION OF NEIGHBOURHOODS
#Select the neighbourhood(s) you want to study
# to see the names of neighbourhoods: unique(all_trees$Neighbourhood)
list.of.neighbourhoods      = list()
list.of.neighbourhoods[[1]] = "NORTH RIVER HEIGHTS"
####################################################################################

#libraries
library(sqldf)
library(igraph)
library(Matrix)
library(R.utils)

## Fuctions

###PROBA_ROOTS
#
#Compute the probability to be infected through the roots using the distance dist
proba_roots = function(dist, minDist, maxDist){
  #p_r proba to infect another tree if dist<min, where min is h(T)
  if(dist<=minDist){
    out = 1
  } else if (dist>=maxDist) {
    out = 0
  } else {
    out = (maxDist-dist)/(maxDist-minDist)
  }
  return(out)
}

###PROBA_ALL_ROOTS
#
#returns a vector that, for each couple of trees, compute the probability to be infected through roots
#Here, the maximum distance at which trees can be connected is given by 3 times the sum of each 
proba_all_roots = function(distances_Neighbourhood){
  n = dim(distances_Neighbourhood)[1]
  vec_proba = mat.or.vec(n,1)
  for (i in 1:n){
    min_d = distances_Neighbourhood$height_i+distances_Neighbourhood$height_j
    max_d = 3*(distances_Neighbourhood$height_i+distances_Neighbourhood$height_j)
    vec_proba[i] = proba_roots(distances_Neighbourhood$dist[i],min_d,max_d)
  }
  return(vec_proba)
}

###SELECT_TREES_NEIGHBOURHOOD
#
#For a given neighbourhood, select elms from the entire dataset all_trees
select_trees_neighbourhood = function(all_trees,Neighbourhood,save_file){
  elms = all_trees[which(all_trees$Neighbourhood == Neighbourhood),]
  saveRDS(elms,file = save_file)
  return(elms)
}

# Set directories
source(sprintf("%s/CODE/set_directories.R", here::here()))

# Put the name of the last tree inventory created by pre_roots_vs_routes.R
all_trees = readRDS(sprintf("%s/tree_inventory_elms_2020-08-26.Rds", DIRS$DATA))

#not needed anymore
elms_idx = grep("American Elm",
                all_trees$Common.Name,
                ignore.case = TRUE)
all_trees = all_trees[elms_idx,]
all_trees = all_trees[which(all_trees$DBH>5),] #remove trees too small
all_trees$Index = 1:dim(all_trees)[1]
distances = readRDS(sprintf("%s/elms_distances_roots.Rds",DIRS$DATA))


#The following loop makes 1) the separation into neighbourhoods for the selected neighbourhoods in list.neighbourhoods, 2) select the right lines in the elm_distances_root that correspond to the trees in the neighbourhood and 3) compute the probabilities that two close trees infect can infect each other through the roots
list_trees = list()
for (Neighbourhood in list.of.neighbourhoods){
  Neighbourhood_norm = gsub(" ", "_", as.character(Neighbourhood), fixed = TRUE)
  Neighbourhood_norm = gsub(".", "", as.character(Neighbourhood_norm), fixed = TRUE)
  Neighbourhood_norm = gsub("'", "", as.character(Neighbourhood_norm), fixed = TRUE)
  #first, save the dataset of trees for each neighbourhood
  save_file = sprintf("%s/Elms_%s.Rds", DIRS$elms, Neighbourhood_norm)
  trees = select_trees_neighbourhood(all_trees, Neighbourhood, save_file)
  list_trees[[Neighbourhood_norm]] = trees
  #second, take the dataframe with distances and select source and destination trees that are in the Neighbourhood (all elms in trees)
  idx_sources = distances$ID_i %in% trees$Tree.ID
  idx_sources = which(idx_sources)
  idx_destinations = distances$ID_j[idx_sources] %in% trees$Tree.ID
  idx_destinations = which(idx_destinations)
  distances_Neighbourhood = distances[idx_sources,][idx_destinations,]
  #third, compute the probabilities for each tree to 
  distances_Neighbourhood$proba    = proba_all_roots(distances_Neighbourhood)
  #and we need to double the dataframe, so that in the simulations, we just look once over sources.
  sub_i = distances_Neighbourhood[,1:8]
  sub_j = distances_Neighbourhood[,9:16]
  sub_k = distances_Neighbourhood[,17:23]
  super_sub = cbind(sub_j,sub_i,sub_k)
  colnames(super_sub) = colnames(distances_Neighbourhood)
  double_distances_Neighbourhood = rbind(distances_Neighbourhood,super_sub)
  saveRDS(distances_Neighbourhood,sprintf("%s/Proba_roots_%s.Rds",DIRS$roots,Neighbourhood_norm))
  #Finally, save the number of trees
  sim_core = dim(trees)[1]
  saveRDS(sim_core,sprintf("%s/sim_core_%s.Rds",DIRS$DATA,Neighbourhood_norm))
}