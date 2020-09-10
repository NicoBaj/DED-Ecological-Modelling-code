### SET MAIN DIRECTORY AND MOVE THERE

# Store all the information in a list
DIRS = list()

DIRS$CODE = "."
DIRS$DATA = "DATA"
DIRS$RESULTS = "RESULTS"

# Where we save some of the processed data (finer than just DATA)
DIRS$neighbourhoods = sprintf("%s/neighbourhoods", DIRS$DATA)
# Some prefixing and suffixing info that will be needed 
DIRS$prefix_data_date = sprintf("%s/data_", DIRS$neighbourhoods)
DIRS$suffix_preproc_dists = "pre_processed_distances"