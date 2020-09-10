# SET_DIRECTORIES.R
#
#
# When sourced, this file should be able to define the location of all the
# needed directories for the various scripts and functions.
# 
# This file is used to produce simulations in the paper: 
# Spread of Dutch Elm Disease in an urban forest
# Nicolas Bajeux, Julien Arino, Stephanie Portet and Richard Westwood
# Ecological Modelling

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