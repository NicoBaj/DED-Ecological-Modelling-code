# PRE_NEIGHBOURHOODS_PROBA_ROOTS.R
#
#
# Load the right dataset for a neighbourhood, compute the probabilities of
# root-based infection.
# 
# This file is used to produce simulations in the paper: 
# Spread of Dutch Elm Disease in an urban forest
# Nicolas Bajeux, Julien Arino, Stephanie Portet and Richard Westwood
# Ecological Modelling

####################################################################################
### SELECTION OF NEIGHBOURHOODS
# Select the neighbourhood(s) you want to study
# To see the names of neighbourhoods: unique(all_trees$Neighbourhood), once you have 
# loaded the variable all_trees (a little below). So you may want to run this script a 
# first time, then explore the contents of that variable.
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
  if(dist<minDist){
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
proba_all_roots = function(distances_nbhd){
  n = dim(distances_nbhd)[1]
  vec_proba = mat.or.vec(n,1)
  for (i in 1:n){
    min_d = distances_nbhd$height_i[i]+distances_nbhd$height_j[i]
    max_d = 3*(distances_nbhd$height_i[i]+distances_nbhd$height_j[i])
    vec_proba[i] = proba_roots(distances_nbhd$dist[i],min_d,max_d)
  }
  return(vec_proba)
}

###SELECT_TREES_NEIGHBOURHOOD
#
#For a given neighbourhood, select elms from the entire dataset all_trees
select_trees_neighbourhood = function(all_trees, nbhd, save_file){
  elms = all_trees[which(all_trees$Neighbourhood == nbhd),]
  saveRDS(elms,file = save_file)
  return(elms)
}

# Set directories
source(sprintf("%s/CODE/set_directories.R", here::here()))

# There can be several versions of the tree inventory file, load the latest
TI_files = list.files(path = DIRS$DATA,
                      pattern = glob2rx("tree_inventory_elms*.Rds"))
selected_TI_file = sort(TI_files, decreasing = TRUE)[1]
# Override if needed by selecting manually one of the files in TI_files, for instance
# selected_TI_file = TI_files[1]
# Get the date, to save distance files with that information
date_TI_file = substr(selected_TI_file, 21, 30)
# Load elms inventory file
all_trees = readRDS(file = sprintf("%s/%s", DIRS$DAT, selected_TI_file))
# Load all distances between these trees
distances = readRDS(sprintf("%s/elms_distances_roots_%s.Rds", DIRS$DATA, date_TI_file))
# Set save directory to include date of data file
DIRS$nbhd_and_date = sprintf("%s%s", DIRS$prefix_data_date, date_TI_file)


# The following loop makes 1) the separation into neighbourhoods for the selected neighbourhoods in 
# list.neighbourhoods, 2) select the right lines in the elm_distances_root that correspond to the 
# trees in the neighbourhood and 3) compute the probabilities that two close trees infect can infect 
# each other through the roots
list_trees = list()
for (nbhd in list.of.neighbourhoods){
  nbhd_norm = gsub(" ", "_", as.character(nbhd), fixed = TRUE)
  nbhd_norm = gsub(".", "", as.character(nbhd_norm), fixed = TRUE)
  nbhd_norm = gsub("'", "", as.character(nbhd_norm), fixed = TRUE)
  # first, save the dataset of trees for each neighbourhood
  save_file = sprintf("%s/elms_%s.Rds", DIRS$nbhd_and_date, nbhd_norm)
  trees = select_trees_neighbourhood(all_trees, nbhd, save_file)
  list_trees[[nbhd_norm]] = trees
  #second, take the dataframe with distances and select source and destination trees that are in the nbhd (all elms in trees)
  idx_sources = distances$ID_i %in% trees$Tree.ID
  idx_sources = which(idx_sources)
  idx_destinations = distances$ID_j[idx_sources] %in% trees$Tree.ID
  idx_destinations = which(idx_destinations)
  distances_nbhd = distances[idx_sources,][idx_destinations,]
  #third, compute the probabilities for each tree to connect to neighbours
  distances_nbhd$proba = proba_all_roots(distances_nbhd)
  #and we need to double the dataframe, so that in the simulations, we just look once over sources.
  sub_i = distances_nbhd[,1:8]
  sub_j = distances_nbhd[,9:16]
  sub_k = distances_nbhd[,17:23]
  super_sub = cbind(sub_j,sub_i,sub_k)
  colnames(super_sub) = colnames(distances_nbhd)
  double_distances_nbhd = rbind(distances_nbhd,super_sub)
  saveRDS(double_distances_nbhd,sprintf("%s/proba_roots_%s.Rds",
                                        DIRS$nbhd_and_date,
                                        nbhd_norm))
  #Finally, save the number of trees
  sim_core = dim(trees)[1]
  saveRDS(sim_core,sprintf("%s/sim_core_%s.Rds", DIRS$nbhd_and_date, nbhd_norm))
}