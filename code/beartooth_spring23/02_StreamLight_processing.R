## Recovery of Stream Productivity following Disturbance
## Originally created: October 12, 2022
## Heili Lowman

#### READ ME ####

# The following set of scripts will walk through the steps necessary to
# prep and send data to Beartooth as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_XX" folders have been ignored using git.ignore.
# Links to the raw data sets are provided in the 01_NWIS_RiverSelection.R file.

# If you are accessing the code via GitHub, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running the 
# code below.

#### Setup ####

lapply(c("tidyverse", "plyr", "cowplot","lubridate",
         "data.table","patchwork", "here"), require, character.only=T)

# Load in datasets created in NWIS_RiverSelection script (Step 1).
NWIS_198sites_subset <- readRDS("data_working/NWIS_198sites_050423.rds")
NWIS_198sitesinfo_subset <- readRDS("data_working/NWIS_198sitesinfo_050423.rds")
NWIS_198sites_Ndays <- readRDS("data_working/NWIS_198sitesNdays_050423.rds")

# Load in list of StreamLight data downloaded from Phil's GitHub repo.
# Note, this provides cumulative daily light, which aligns with the 
# magnitude of light values used in the original Appling modeling efforts.
# https://github.com/streampulse/metabolism_synthesis/tree/master/output_data
# "Stream_PAR_sum": Daily sum of PAR estimated at the stream surface (mol m^-2^ d^-1^)
rawSL <- readRDS("data_raw/metabolism_synthesis_output_data/lotic_standardized_full.rds")

#### Trim, evaluate NAs, and trim again ####

# Savoy SL data
SL_Sav_df <- plyr::ldply(rawSL, data.frame)

# trim down to necessary columns
SL_Sav_df <- SL_Sav_df %>%
  select(Site_ID, Date, DOY, Year, Stream_PAR_sum) 

# and rename columns to match other datasets
SL_Sav_df <- SL_Sav_df %>%
  dplyr::rename(site_name = Site_ID,
         PAR_surface = Stream_PAR_sum)

# Remove NA and NaN values, because these will cause the model to break.
SL_Sav_clean <- SL_Sav_df %>%
  filter(!is.na(PAR_surface)) %>% # removes NAs
  filter(!is.nan(PAR_surface)) # removes NaNs
# 804,488 observations remaining

SL_Sav_trim <- SL_Sav_clean %>%
  filter(site_name %in% NWIS_198sitesinfo_subset$site_name)
# 449,150 observations remaining

# And create a list form as well.
SL_trim_list <- split(SL_Sav_trim, SL_Sav_trim$site_name)

## Save both formats, df and list.
saveRDS(SL_Sav_trim, "data_working/StreamLight_daily_df_182sites_050423.rds")
saveRDS(SL_trim_list, "data_working/StreamLight_daily_list_182sites_050423.rds")

# End of script.
