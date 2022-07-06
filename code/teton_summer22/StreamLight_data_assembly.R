## Remaining StreamLight data source
## July 5, 2022
## Heili Lowman

# Datasets downloaded from : https://github.com/streampulse/metabolism_synthesis/tree/master/output_data

## Load packages
lapply(c("tidyverse", "lubridate", "rstan",
         "rlist", "pipeR", "here"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################

# Load in list of sites for which StreamLight data is still missing.
siteswoSL <- read_csv("data_working/comids_100sites.csv")

# Load in list of StreamLight data downloaded from Phil's GitHub repo.
rawSL <- readRDS("data_raw/metabolism_synthesis_output_data/lotic_standardized_full.rds")

# Create list to filter by.
siteswoSL <- siteswoSL %>%
  mutate(nwis_ID = paste('nwis_', nwis_id, sep = ""))

mylist <- siteswoSL$nwis_ID

# Filter dataset by list of sites needing StreamLight data.

# Create custom function
myfilter <- function(x){
  
  # filter by list of sites created above
  df <- x %>%
    filter(Site_ID %in% mylist)
  
  # return new dataset (if matched)
  df
  
}

filterSL <- lapply(rawSL, myfilter) # only 7 sites matched :(

# End of script.


