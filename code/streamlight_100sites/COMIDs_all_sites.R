## 207 rivers data source
## StreamLight Data Acquisition
## June 10, 2022
## Heili Lowman

# Pulling out COMIDs for azimuth calculation at all sites.

# Load packages
lapply(c("tidyverse", "here",
         "StreamLightUtils", "StreamLight"), require, character.only=T)

# Importing site list displaying which sites have StreamLight data available.
sitelist <- read_csv("data_working/streamlight_sites.csv")

# Importing site info list with lat/lon information
site_info <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")

# Need to create "A table with Site Name, COMID, Lat, and Lon, and Long Name"
# so, filter for sites that don't have Stream Light data already available
sites_need_SL <- sitelist %>%
  dplyr::filter(SL == "NO")

# then, combine the two datasets
my_site_comids <- left_join(sites_need_SL, site_info, by = c("sites" = "site_name")) %>%
  dplyr::select(-SL) # remove stream light indicator column

comids_trim <- my_site_comids %>%
  dplyr::select(nwis_id, nhdplus_id, lat, lon, long_name)

write_csv(comids_trim, "data_working/comids_100sites.csv")

# End of script.
