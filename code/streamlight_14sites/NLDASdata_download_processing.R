## 207 rivers data source
## StreamLight Data Acquisition
## June 2, 2022
## Heili Lowman

# Working through first step of Stream Light data download workflow described at:
# https://psavoy.github.io/StreamLight/articles/1%20Download%20and%20process%20NLDAS.html

# "This article describes downloading and processing total incoming shortwave radiation (W m-2) using the StreamLightUtils package. Total incoming shortwave radiation (W m-2) is provided by the National Land Data Assimilation System (NLDAS) at hourly timesteps at 0.125 degree spatial resolution."

## Load packages
lapply(c("tidyverse", "here",
         "StreamLightUtils", "StreamLight"), require, character.only=T)

## Load datasets.

# Importing site list displaying which sites have StreamLight data available.
sitelist <- read_csv("data_working/streamlight_sites.csv")

# Importing site info list with lat/lon information
site_info <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")

# Need to create "A table with Site_ID, Lat, and Lon, and startDate"
# so, filter for sites that don't have Stream Light data already available
sites_need_SL <- sitelist %>%
  dplyr::filter(SL == "NO")

# and filter for necessary geographic information
site_lls <- site_info %>%
  dplyr::select(site_name, lat, lon)

# then, combine the two datasets
my_site_locs <- left_join(sites_need_SL, site_lls, by = c("sites" = "site_name")) %>%
  dplyr::select(-SL) %>%
  dplyr::rename("Site_ID" = "sites",
         "Lat" = "lat",
         "Lon" = "lon") %>%
  mutate(startDate = "2007-01-01")

# The following has been edited from the original code in the documentation:

#Set the download location (add your own directory)
working_dir <- "C:/Users/hlowman/Documents/ContinentalRiverBiomass/data_raw/StreamLight_NLDAS"

#Read in a table with initial site information
sites <- my_site_locs

#Download NLDAS data at 100 sites for which data was not published
NLDAS_DL_bulk(
  save_dir = "C:/Users/hlowman/Documents/ContinentalRiverBiomass/data_raw/StreamLight_NLDAS",
  site_locs = sites
) #WHAT!?! These files appear and then immediately dissappear...

# This isn't working, so trying it at a single site to troubleshoot
NLDAS_DL("C:/Users/hlowman/Documents/ContinentalRiverBiomass/data_raw/StreamLight_NLDAS", "nwis_01124000", 42.02232, -71.95563, "2007-01-01")
# Ok, so this works.

# Instead going to build a loop to loop over sites and see if this sticks.
for (i in 1:100) {
  NLDAS_DL("C:/Users/hlowman/Documents/ContinentalRiverBiomass/data_raw/StreamLight_NLDAS", my_site_locs$Site_ID[i], my_site_locs$Lat[i], my_site_locs$Lon[i], "2007-01-01")
}

# Alright, this seems to have worked for 98 sites.

#List of successfully downloaded sites
NLDAS_list <- stringr::str_sub(list.files("C:/Users/hlowman/Documents/ContinentalRiverBiomass/data_raw/StreamLight_NLDAS"), 1, -11)

#Processing the downloaded NLDAS data
NLDAS_processed <- StreamLightUtils::NLDAS_proc(read_dir = "C:/Users/hlowman/Documents/ContinentalRiverBiomass/data_raw/StreamLight_NLDAS", NLDAS_list)

saveRDS(NLDAS_processed, "C:/Users/hlowman/Documents/ContinentalRiverBiomass/data_raw/NLDAS_processed_072822.rds")

# End of NLDAS data download and processing script.
