## Resilience of Stream Productivity to Disturbance
## November 30, 2022
## Heili Lowman


#### Setup ####

# Load necessary packages.
lapply(c("tidyverse", "lubridate", "data.table",
         "nhdplusTools", "sf", "mapview", "here"), require, character.only=T)

#### Data ####

# Trying out get_huc() function.
testhuc <- get_huc(id = 06144833, buffer = 0.5, type = "huc10")

# Spits back this:
# Spherical geometry (s2) switched off
# Spherical geometry (s2) switched on
# Warning message:
#   No huc10 features found 

# Ok, so the COMID approach doesn't work.

# Trying Mike V.'s workflow

location <- data.frame(x = 42.02232, y = -71.95563) %>%
  st_as_sf(coords = c('y', 'x'), crs = 4326)

huc12 = nhdplusTools::get_huc(location, type = "huc12") # only works one site at a time

mapview(huc12) + mapview(location) # map-truth site + watershed

huc10_id = substr(huc12$huc12, 1, 10) # OMG WUT

# Ok, so load in rmax data to get the 159 sites for which I need HUC10 data.
dat_rmax <- readRDS("data_working/rmax_filtered_159sites_113022.rds")

# And trim down to sites and their lat/lon.
dat_geo <- dat_rmax %>%
  dplyr::select(site_name, Lat_WGS84, Lon_WGS84)

# Make new dataset with sf geometries.
dat_locations <- dat_geo %>%
  rename(x = Lat_WGS84, y = Lon_WGS84) %>%
  st_as_sf(coords = c('y', 'x'), crs = 4326)

huc12 <- vector("list", 159) # create an empty list

for(i in 1:length(dat_locations$site_name)){
  
  huc12[[i]] <- nhdplusTools::get_huc(dat_locations$geometry[i], type = "huc12")
  
}

# Testing a few maps to see if this worked properly.

mapview(huc12[10]) + mapview(dat_locations$geometry[10])

mapview(huc12[150]) + mapview(dat_locations$geometry[150])

mapview(huc12[90]) + mapview(dat_locations$geometry[90])

# Well shoot, this is amazing.

# Make the list into a dataframe.
huc12_df <- do.call(rbind, huc12)

# And add huc10 to it.
huc12_df$huc10_id <- substr(huc12_df$huc12, 1, 10)

# Quick bar plot to visualize how many sites are in the same watershed.
ggplot(huc12_df, aes(x = huc10_id)) + geom_bar()

# So, roughly one third of sites share a watershed with other sites.

# And export this data for use in the linear-mixed effects modeling.
saveRDS(huc12_df, file = "data_working/HUC12_159sites_120722.rds")

# End of script.
