## 207 rivers data source
## StreamLight Data Acquisition
## July 28, 2022
## Heili Lowman

# The following script will pull together an azimuth dataset for the 14 remaining sites
# for which I need to calculate light data.

# Load packages
lapply(c("tidyverse", "here",
         "StreamLightUtils", "StreamLight"), require, character.only=T)

# Create a list of the sites of interest.
my14sites <- c("nwis_023362095","nwis_01656903","nwis_04174500",
"nwis_385520094420000","nwis_06893300","nwis_04168400","nwis_03098600",
"nwis_04167150","nwis_02336120","nwis_02336240","nwis_02336152",
"nwis_08057000","nwis_04199500","nwis_040871488")

# Read in site info to get locations.
NWIS_207sitesinfo_subset <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")

# Filter down to only the sites of interest.
NWIS_14sitesinfo <- NWIS_207sitesinfo_subset %>%
  dplyr::filter(site_name %in% my14sites)

# And trim the dataset down some.
NWIS_14sitesinfo <- NWIS_14sitesinfo[,1:6]

# Using Google Maps, and the instructions for deriving azimuth at each
# site, I determined azimuth values manually for the remaining 14 sites.
# https://psavoy.github.io/StreamLight/articles/3%20Using%20stream_light.html

# Workflow: type lat + lon into Google Maps search bar, and then zoom in 3x
# prior to determining approximate azimuth value.

my14azimuths <- c(45, 315, 80, 45, 345, 315, 350, 65, 55, 45, 35, 95, 345, 5)

# Add this info to the dataset and save it out.
NWIS_14sitesinfo <- NWIS_14sitesinfo %>%
  mutate(azimuth = my14azimuths)

saveRDS(NWIS_14sitesinfo, "C:/Users/hlowman/Documents/ContinentalRiverBiomass/data_raw/azimuths_14sites_072822.rds")

# End of script.
