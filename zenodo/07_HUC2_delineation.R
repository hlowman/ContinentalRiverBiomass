## Recovery of Stream Productivity following Disturbance
## Originally created: November 30, 2022
## Heili Lowman

#### READ ME ####

# The following set of scripts will walk through the steps necessary to
# prep and send data to the University of Wyoming's Beartooth High Performance
# Computing Cluster as well as process the model outputs.

# Please note, the "data_XX" folders have been ignored using git.ignore.
# Links to the raw data sets are provided in the 01_NWIS_RiverSelection.R file
# and 02_StreamLight_processing files.

# If you are accessing the code via GitHub or Zenodo, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running the 
# code below.

# This script adds information regarding the HUC2 identifier for each site.

#### Setup ####

# Load necessary packages.
lapply(c("tidyverse", "lubridate", "data.table", "mapview",
         "nhdplusTools", "sf", "here"), require, character.only=T)

# Note, the mapview package can be VERY finicky.

#### Data ####

# Using Mike V.'s workflow

location <- data.frame(x = 42.02232, y = -71.95563) %>%
  st_as_sf(coords = c('y', 'x'), crs = 4326)

huc12 = nhdplusTools::get_huc(location, type = "huc12") # only works one site at a time

mapview(huc12) + mapview(location) # map-truth site + watershed

huc10_id = substr(huc12$huc12, 1, 10) # yay!

# Ok, so load in NWIS data to get the sites for which I need HUC10 data.
# Can't use the HUCs in this dataset because they are HUC8.
dat <- readRDS("data_working/NWIS_Info_181riv_050923.rds")

# And trim down to sites and their lat/lon.
dat_df <- map_df(dat, ~as.data.frame(.x) %>%
                   mutate(map_scale_fc = as.character(map_scale_fc)), .id="site_name") 

dat_geo <- dat_df %>%
  dplyr::select(site_no, dec_lat_va, dec_long_va)

# Make new dataset with sf geometries.
dat_locations <- dat_geo %>%
  rename(x = dec_lat_va, y = dec_long_va) %>%
  st_as_sf(coords = c('y', 'x'), crs = 4269)

huc12 <- vector("list", 181) # create an empty list

for(i in 1:length(dat_locations$site_no)){
  
  huc12[[i]] <- nhdplusTools::get_huc(dat_locations$geometry[i], type = "huc12")
  
}

# Testing a few maps to see if this worked properly.

mapview(huc12[30]) + mapview(dat_locations$geometry[30])

mapview(huc12[150]) + mapview(dat_locations$geometry[150])

mapview(huc12[90]) + mapview(dat_locations$geometry[90])

# Well shoot, this is amazing.

# Make the list into a dataframe.
huc12_df <- do.call(rbind, huc12)

# And add huc10 to it.
huc12_df$huc10_id <- substr(huc12_df$huc12, 1, 10)

# and add site names to it.
huc12_df$site_no <- dat_locations$site_no

# Quick bar plot to visualize how many sites are in the same watershed.
ggplot(huc12_df, aes(x = huc10_id)) + geom_bar()

# So, roughly one third of sites share a watershed with other sites.

#### Additional HUC groupings ####

# Since the HUC10 level had too few sites per grouping (mean = 1), I'm
# deriving the remaining data from the above dataset.

#huc12_df <- readRDS("data_working/HUC12_159sites_120722.rds")

huc12_df$huc8_id <- substr(huc12_df$huc12, 1, 8)
huc12_df$huc6_id <- substr(huc12_df$huc12, 1, 6)
huc12_df$huc4_id <- substr(huc12_df$huc12, 1, 4)
huc12_df$huc2_id <- substr(huc12_df$huc12, 1, 2)

# And see how many groups there are for each (looking to hit ~15/group)
ggplot(huc12_df, aes(x = huc8_id)) + geom_bar() # still mostly 1/group
ggplot(huc12_df, aes(x = huc6_id)) + geom_bar() # a bit better but ~2/group
ggplot(huc12_df, aes(x = huc4_id)) + geom_bar() # not a single group at 15 yet
ggplot(huc12_df, aes(x = huc2_id)) + geom_bar() # ok now we might be talking

# So, after all that work, I'm simply going to append a HUC2 column to the
# existing NWIS info dataset.

dat_df <- dat_df %>%
  mutate(huc_2 = substr(huc_cd, start = 1, stop = 2))

# Export data.
saveRDS(dat_df, file = "data_working/NWIS_Info_181riv_HUC2_df_050923.rds")

# End of script.
