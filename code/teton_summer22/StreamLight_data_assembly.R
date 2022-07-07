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

# Ok, so looks like I'll need to work back through the StreamLight workflow.

# Just kidding! Sometimes an extra 0 is added to site names, so workflow below should take care of that and re-match the dataframes.

# Using Phil's code to pull the sites I need
sites_100 <- siteswoSL

#Merge together the site information
info_merged <- merge(
  sites_100[, c("nwis_id", "nhdplus_id", "lat", "lon", "long_name")], 
  lotic_site_info_full, 
  by.x = "long_name",
  by.y = "Name"
)

#Subset the list just for the sites of interest
standard_ized_subset <- lotic_standardized_full[info_merged[, "Site_ID"]]

#Check to see if there are no predictions of light at the stream surface
light_check <- do.call(
  rbind,
  lapply(
    standard_ized_subset,
    FUN = function(site){all(is.na(site[, "Stream_PAR_sum"]))}
  )  
) %>%
  as.data.frame()

light_check$Site_ID <- rownames(light_check)
names(light_check)[1] <- "no_light"

#Get string of sites that do have light estimates
has_light <- light_check %>%
  filter(no_light == FALSE) %>%
  pull(Site_ID)

# End of Phil's code.

# Subset the list just for the sites with stream light data
final_subset <- standard_ized_subset[has_light]

# Export data
saveRDS(final_subset, "data_working/StreamLight_daily_85site_list.rds")

# End of script.


