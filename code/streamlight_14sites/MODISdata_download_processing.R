## 207 rivers data source
## StreamLight Data Acquisition
## June 2, 2022
## Heili Lowman

# Working through second step of Stream Light data download workflow described at:
# https://psavoy.github.io/StreamLight/articles/2%20Download%20and%20process%20MODIS%20LAI.html

# "There are a variety of sources for LAI data, but for convenience a function is included in StreamLightUtils to process two MODIS LAI products: 1.) MCD15A3H.006 which has 4 day temporal resolution and a pixel size of 500m, and 2.) MCD15A2H.006 which has 8 day temporal resolution and a pixel size of 500m."

# Load packages
lapply(c("tidyverse", "here",
         "StreamLightUtils", "StreamLight"), require, character.only=T)

# Load datasets & process similarly to NLDAS script.

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
sites <- my_site_locs

#Make a table for the MODIS request 
request_sites <- sites[, c("Site_ID", "Lat", "Lon")]

#Export your sites as a .csv for the AppEEARS request  
write.table(
  request_sites, 
  paste0("data_working/PM_sites.csv"), 
  sep = ",", 
  row.names = FALSE,
  quote = FALSE, 
  col.names = FALSE
)

# Request for data mades through AppEEARS website at 2:30pm 6/2/2022.

# SECOND REQUEST: 11:50am 7/28/2022.
# I've re-run the request because the previous output wasn't cooperating
# when I noticed an error in how I'd saved things out and tried to re-do it.

# Create a list of the sites of interest.
my14sites <- c("nwis_023362095","nwis_01656903","nwis_04174500",
               "nwis_385520094420000","nwis_06893300","nwis_04168400",
               "nwis_03098600","nwis_04167150","nwis_02336120",
               "nwis_02336240","nwis_02336152",
               "nwis_08057000","nwis_04199500","nwis_040871488")

# New trimmed dataset for second request.
request_sites2 <- request_sites %>%
  dplyr::filter(Site_ID %in% my14sites) %>%
  dplyr::select(Site_ID, Lat, Lon)

#Export your sites as a .csv for the AppEEARS request  
write.table(
  request_sites2, 
  paste0("data_working/PM_sites14.csv"), 
  sep = ",", 
  row.names = FALSE,
  quote = FALSE, 
  col.names = FALSE
)

# When I went back to re-do this, I had to delete the former PM_site folder, 
# reset the working directory using the code below, and re-run the following.
# setwd("C:/Users/hlowman/Documents/ContinentalRiverBiomass/data_raw/StreamLight_MODIS")

# Also, the AppEEARS website drops the "_" in all the sitenames, so I
# need to create a new column that does that as well, otherwise no data
# is pulled out of the .csv file.
request_sites2 <- request_sites2 %>%
  mutate(SiteID = gsub("_","", Site_ID))

# Unpack MODIS data received.
MOD_unpack <- AppEEARS_unpack_QC(
  zip_file = "PM_sites14.zip", 
  zip_dir = "C:/Users/hlowman/Documents/ContinentalRiverBiomass/data_raw/StreamLight_MODIS", 
  request_sites2[, "SiteID"]
)

# Process and plot MODIS data received.
MOD_processed <- AppEEARS_proc(
  unpacked_LAI = MOD_unpack,  
  fit_method = "Gu", 
  plot = TRUE
)

# Received the following error:
# Error in ts[, "FparLai_QC_SCF_QC"] : incorrect number of dimensions

# saveRDS(MOD_processed, "C:/Users/hlowman/Documents/ContinentalRiverBiomass/data_raw/MODIS_processed_072822.rds")

# End of MODIS data download and processing script.
