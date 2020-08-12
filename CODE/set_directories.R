### SET MAIN DIRECTORY AND MOVE THERE

# Name of the nodes on gauss
nodes = c("node100","node101","node102","node103")

# DATA AND OUTPUT
# If we are working on the cluster, use absolute path, otherwise use relative to the user
if (Sys.info()["nodename"] %in% nodes) { 
  TOP_DIR_DATA = "/storage/home/bajeuxni/DED-Ecological-Modelling-code/DATA"
  TOP_DIR_RES  = "/storage/home/bajeuxni/DED-Ecological-Modelling-code/RESULTS"
  # TOP_DIR_DATA_OUTPUT = "/storage/var/groups/mathbio/DED_DATA_OUTPUT"
  TOP_DIR_CODE = "/storage/home/bajeuxni/DED-Ecological-Modelling-code/CODE"
} else if (Sys.info()["nodename"] == "DESKTOP-OC4RBTE") {
  TOP_DIR_DATA_OUTPUT = "C:/Users/nicol/Dropbox/Postdoc_Canada/winnipeg_trees/Shiny/DED/code"
} else {#for the shiny
  TOP_DIR_DATA_OUTPUT = "./code"
}

# Set working directory at top of data and output
setwd(TOP_DIR_DATA_OUTPUT)
