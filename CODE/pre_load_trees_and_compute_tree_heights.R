# PRE_LOAD_TREES_AND_COMPUTE_TREE_HEIGHTS.R
#
#
# Load tree data, select elms and add tree heights.
# 
# This file is used to produce simulations in the paper: 
# Nicolas Bajeux, Julien Arino, Stephanie Portet and Richard Westwood

# Set directories
source(sprintf("%s/CODE/set_directories.R", here::here()))

# Tree height as a function of diametre at breast height
th <- function(dbh) {
  out <- round(4.208*log(dbh)-1.707,2)
  return(out)
}

# Do we want the latest tree data?
REFRESH_ELM_DATA = TRUE

if (REFRESH_ELM_DATA) {
  Tree_Inventory <- read.csv("https://data.winnipeg.ca/api/views/hfwk-jp4h/rows.csv?accessType=DOWNLOAD",
                             stringsAsFactors = FALSE)
} else {
  Tree_Inventory <- read.csv(sprintf("%s/Tree_Inventory.csv",DIRS$DATA),
                             stringsAsFactors = FALSE)
}

# Select American Elms
elms = Tree_Inventory[grep("American Elm", 
                           Tree_Inventory$Common.Name,
                           ignore.case = TRUE),]
# We want to change the column name for DBH (to this)
names_cols_elms = names(elms)
idx_dbh = which(names_cols_elms == "Diameter.at.Breast.Height")
colnames(elms)[idx_dbh] = "DBH"

# Get rid of tiny elms (DBH<=5cm)
elms = elms[which(elms$DBH > 5),]

# Add a column with tree heights based on formula in paper
elms = cbind(elms[,1:idx_dbh],
             th(elms$DBH),
             elms[,(idx_dbh+1):length(names_cols_elms)])
colnames(elms)[(idx_dbh+1)] = "TreeHeight"

# We also need to break up location as lat,lon
elms$lat = lapply(strsplit(elms$Location, ","), function(x) as.numeric(substr(x[1], 2, nchar(x[1]))))
elms$lon = lapply(strsplit(elms$Location, ","), function(x) as.numeric(substr(x[2], 1, (nchar(x[2])-1))))

# Save the resulting file, suffixing with today's date
saveRDS(elms, file = sprintf("%s/tree_inventory_elms_%s.Rds", DIRS$DATA, Sys.Date()))
