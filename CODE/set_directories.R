### SET MAIN DIRECTORY AND MOVE THERE

# Store all the information in a list
DIRS = list()

# Provided this file is sourced, this should set to the top of the Github repo
DIRS$prefix_github = here::here()

DIRS$CODE = sprintf("%s/CODE", DIRS$prefix_github)
DIRS$DATA = sprintf("%s/DATA", DIRS$prefix_github)
DIRS$RESULTS = sprintf("%s/RESULTS", DIRS$prefix_github)

# Where we save some of the processed data (finer than just DATA)
DIRS$neighbourhoods = sprintf("%s/neighbourhoods", DIRS$DATA)
# Some prefixing and suffixing info that will be needed 
DIRS$prefix_data_date = sprintf("%s/data_", DIRS$neighbourhoods)
DIRS$suffix_preproc_dists = "pre_processed_distances"