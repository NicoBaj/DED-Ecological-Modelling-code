# Compute network for roots, taking into account the road network and other obstacles.

library(osmdata)
library(sf)
library(rprojroot)

# Set top directory for code
TOP_DIR_CODE = here::here()
if (length(TOP_DIR_CODE) == 0) {
  TOP_DIR_CODE = rprojroot::find_root(is_git_root,".")
}
# Source code that sets directories. Always run.
source(sprintf("%s/set_directories.R", TOP_DIR_CODE))
# Source some useful functions
source(sprintf("%s/useful_functions.R", TOP_DIR_CODE))

# If you want to refresh the OSM data, set this to TRUE. Otherwise, pre-saved data is used
REFRESH_OSM_DATA = FALSE

if (REFRESH_OSM_DATA) {
  # Get exact bounding polygon for Winnipeg
  bb_poly = getbb(place_name = "winnipeg", format_out = "polygon")
  # Road types to download from OSM
  road_types = c("motorway","trunk","primary","secondary","tertiary","residential","unclassified")
  # List to store all different road types in
  ROADS = list()
  # Get roads of the given type
  for (rt in road_types) {
    ROADS[[rt]] <- opq(bbox = bb_poly) %>%
      add_osm_feature(key = 'highway', value = rt) %>%
      osmdata_sf () %>%
      trim_osmdata (bb_poly)
  }
  # Variable to store union of all roads
  roads = c(ROADS[["motorway"]],
            ROADS[["trunk"]],
            ROADS[["primary"]],
            ROADS[["secondary"]],
            ROADS[["tertiary"]],
            ROADS[["residential"]],
            ROADS[["unclassified"]])
  saveRDS(roads,sprintf("%s/Winnipeg_roads.Rds",DIRS$DATA))
  #same thing with the rail
  rail <- opq(bbox = bb_poly) %>%
    add_osm_feature(key = 'railway', value = "rail") %>%
    osmdata_sf () %>%
    trim_osmdata (bb_poly)
  saveRDS(rivers,sprintf("%s/Winnipeg_rail.Rds",DIRS$DATA))
  #same thing with the rivers
  rivers <- opq(bbox = bb_poly) %>%
    add_osm_feature(key = 'waterway', value = "river") %>%
    osmdata_sf () %>%
    trim_osmdata (bb_poly)
  saveRDS(rivers,sprintf("%s/Winnipeg_rivers.Rds",DIRS$DATA))
  #same thing with the parking lots
  parkings <- opq(bbox = bb_poly) %>%
    add_osm_feature(key = 'amenity', value = "parking") %>%
    osmdata_sf () %>%
    trim_osmdata (bb_poly)
  saveRDS(parkings,sprintf("%s/Winnipeg_parkings.Rds",DIRS$DATA))
  # All sources of root cuts
  all_root_cutters = c(roads, rail, rivers, parkings)
  saveRDS(all_root_cutters,sprintf("%s/Winnipeg_all_root_cutters.Rds",DIRS$DATA))
} else {
  roads = readRDS(sprintf("%s/Winnipeg_roads.Rds",DIRS$DATA))
  rail = readRDS(sprintf("%s/Winnipeg_rail.Rds",DIRS$DATA))
  rivers = readRDS(sprintf("%s/Winnipeg_rivers.Rds",DIRS$DATA))
  parkings = readRDS(sprintf("%s/Winnipeg_parkings.Rds",DIRS$DATA))
  all_root_cutters = readRDS(sprintf("%s/Winnipeg_all_root_cutters.Rds",DIRS$DATA))
}

# There can be several versions of the tree inventory file, load the latest
TI_files = list.files(path = DIRS$DATA,
                      pattern = glob2rx("Tree_Inventory_Elms_*.csv"))
latest_TI_file = sort(TI_files, decreasing = TRUE)[1]
date_of_file = substr(latest_TI_file,21,30)

# Read elms csv file (could also read the RDS..)
elms <- read.csv(sprintf("%s/%s", DIRS$DATA, latest_TI_file),
                 stringsAsFactors = FALSE)

# Compute distances and select the ones matching the criterion. 
# Work with X,Y (which are in metres), rather than lon,lat
elms_xy = cbind(elms$X, elms$Y)
# CAREFUL: The next call returns a large object (>10GB). Only run on a machine with enough memory.
D_dist = dist(elms_xy)
# CAREFUL AGAIN: the next call returns a >20GB object and further requires close to 91GB RAM to work.
D_mat = as.matrix(D_dist)
# Clean up and do garbage collection (force return of memory to the system)
rm(D_dist)
gc()

# Take a very conservative upper bound for root system extent: 6 times the maximum height
elms_max_height = 6*max(elms$TreeHeight)
idx_D_mat = which(D_mat > elms_max_height)
D_mat[idx_D_mat] = 0
# Clean up and garbage collection again..
rm(idx_D_mat)
gc()

saveRDS(D_mat, file = sprintf("%s/matrix_tmp.Rds", DIR_DATA_PROCESSED_GLOBAL))
# Keep only pairs with nonzero distance (i.e., <= 6*elms_max_height).
indices = which(D_mat !=0, arr.ind = TRUE)
# Also, only keep one of the edges, not both directions.
indices = indices[which(indices[,"row"] > indices[,"col"]),]

# Make data frame
DISTS = data.frame(idx_i = indices[,1],
                   ID_i = elms$Tree.ID[indices[,1]],
                   height_i = elms$TreeHeight[indices[,1]],
                   x_i = elms$X[indices[,1]],
                   y_i = elms$Y[indices[,1]],
                   lat_i = elms$lat[indices[,1]],
                   lon_i = elms$lon[indices[,1]],
                   ngbhd_i = elms$Neighbourhood[indices[,1]],
                   idx_j = indices[,2],
                   ID_j = elms$Tree.ID[indices[,2]],
                   height_j = elms$TreeHeight[indices[,2]],
                   x_j = elms$X[indices[,2]],
                   y_j = elms$Y[indices[,2]],
                   lat_j = elms$lat[indices[,2]],
                   lon_j = elms$lon[indices[,2]],
                   ngbhd_j = elms$Neighbourhood[indices[,2]],
                   dist = D_mat[indices])
# Clean up and collect garbage..
rm(D_mat)
gc()

# Save as both csv and RDS
write.csv(DISTS, file = sprintf("%s/elms_distances_roots_%s.csv",DIRS$DATA, date_of_file))
saveRDS(DISTS, file = sprintf("%s/elms_distances_roots_%s.Rds",DIRS$DATA, date_of_file))


# The locations of the origins of the pairs
tree_locs_orig = cbind(DISTS$lon_i, DISTS$lat_i)
# The locations of the destinations of the pairs
tree_locs_dest = cbind(DISTS$lon_j, DISTS$lat_j)

tree_pairs = do.call(sf::st_sfc,
                     lapply(
                       1:nrow(tree_locs_orig),
                       function(i){
                         sf::st_linestring(
                           matrix(
                             c(tree_locs_orig[i,],
                               tree_locs_dest[i,]), 
                             ncol=2,
                             byrow=TRUE)
                         )
                       }
                     )
)

if (FALSE) {
  pdf(file = sprintf("%s/OUTPUT/PDF/pairs_preproc_%s.pdf", TOP_DIR_DED, date_of_file),
      width = 50, height = 50)
  plot(tree_pairs)
  dev.off()
}


# Set both crs to be the same
st_crs(tree_pairs) = st_crs(all_root_cutters$osm_lines$geometry)

tictoc::tic()
tree_pairs_all_root_cutters = sf::st_intersects(x = all_root_cutters$osm_lines$geometry,
                                                y = tree_pairs)
tictoc::toc()

tree_pairs_all_root_cutters_intersect = c()
for (i in 1:length(tree_pairs_all_root_cutters)) {
  if (length(tree_pairs_all_root_cutters[[i]])>0) {
    tree_pairs_all_root_cutters_intersect = c(tree_pairs_all_root_cutters_intersect,
                                              tree_pairs_all_root_cutters[[i]])
  }
}
tree_pairs_all_root_cutters_intersect = sort(tree_pairs_all_root_cutters_intersect)
to_keep = 1:dim(tree_locs_orig)[1]
to_keep = setdiff(to_keep,tree_pairs_all_root_cutters_intersect)

if(FALSE){
  pdf(file = sprintf("%s/OUTPUT/PDF/pairs_postproc_%s.pdf", TOP_DIR_DED, date_of_file),
      width = 50, height = 50)
  plot(tree_pairs[to_keep_roads_rivers])
  dev.off() 
}

# Finally, indicate which distance classes the tree combinations satisfy
h_tmp = mat.or.vec(nr = dim(DISTS)[1], nc =5)
for (h in 1:5) {
  s_tmp = (DISTS$height_i+DISTS$height_j)*h
  h_tmp[which(DISTS$dist <= s_tmp),h] = 1
}
colnames(h_tmp) = sprintf("h%d",1:5)
DISTS = cbind(DISTS, h_tmp)

# Keep only the edges not intersected by a road or a river
DISTS = DISTS[to_keep,]

# Save as both csv and RDS
write.csv(DISTS, file = sprintf("%s/elms_distances_roots_%s.csv",DIR_DATA_PROCESSED_GLOBAL,date_of_file))
saveRDS(DISTS, file = sprintf("%s/elms_distances_roots_%s.Rds",DIR_DATA_PROCESSED_GLOBAL,date_of_file))
