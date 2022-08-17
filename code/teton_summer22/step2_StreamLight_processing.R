## Resilience of Stream Productivity to Disturbance
## July 25, 2022
## Heili Lowman

# The following set of scripts will walk through the steps necessary to
# prep and send data to Teton as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_raw" and "data_working" folders have been ignored
# using git.ignore, so links to the raw data sets are provided in the step1
# file. If you are accessing the code via GitHub, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running the 
# code.

############################
## Setup
############################

lapply(c("tidyverse", "plyr", "cowplot","lubridate",
         "data.table","patchwork", "here"), require, character.only=T)

# Load in datasets created in NWIS_RiverSelection script.
NWIS_207sites_subset <- readRDS("data_working/NWIS_207sites_subset2022.rds")
NWIS_207sitesinfo_subset <- readRDS("data_working/NWIS_207sitesinfo_subset2022.rds")
NWIS_207sites_Ndays <- readRDS("data_working/NWIS_207sites_Ndays2022.rds")

# Load in dataset created in StreamLight_data_assembly script.
SL_85sites <- readRDS("data_working/StreamLight_daily_85site_list.rds")

############################
## Join Allison + Phil's datasets
############################

## Sites of interest that have StreamLight data
sites <- NWIS_207sitesinfo_subset$site_name
sites_files <- rep(NA, length(sites))
for(i in 1:length(sites)){
  sites_files[i] <- paste(sites[i],"_input_output.txt", sep="")
}

# need to identify the intersection of files with available light data
filenames <- intersect(sites_files, list.files(here("data_raw/individual_files")))

## Import StreamLight
SL <- plyr::ldply(filenames, function(filename) {
  d <- read.table(here("data_raw", "individual_files", filename), header = T, sep = "\t")
  d$file <- filename
  return(d)
})

## take the mean daily incoming PAR at the surface
SL_split <- split(SL, SL$file)

meandaily_PAR <- function(y){
  df <- y %>%
  group_by(jday) %>%
  summarize_at(.vars = c("DOY","Year","PAR_surface","PAR_turb"), .funs = mean, na.rm = TRUE)

  df$origin <- as.Date(paste0(df$Year, "-01-01"),tz = "UTC") - days(1)
  df$Date <- as.Date(df$DOY, origin = df$origin, tz = "UTC")

  return(df)
}

SL_daily <- lapply(SL_split, function(x) meandaily_PAR(x))

###########################
## Join and evaluate NAs
###########################
# Appling SL data
SL_App_df <- plyr::ldply(SL_daily, data.frame)

## add site_name
SL_App_df$site_name <- substr(SL_App_df$.id, 1, nchar(SL_App_df$.id)-17)

# trim down to necessary columns
SL_App_df <- SL_App_df %>%
  select(site_name, Date, DOY, Year, PAR_surface)


# Savoy SL data
SL_Sav_df <- plyr::ldply(SL_85sites, data.frame)

# trim down to necessary columns
SL_Sav_df <- SL_Sav_df %>%
  select(Site_ID, Date, DOY, Year, Stream_PAR_sum) 

# and rename columns to match the trimmed Appling dataset
SL_Sav_df <- SL_Sav_df %>%
  dplyr::rename(site_name = Site_ID,
         PAR_surface = Stream_PAR_sum)

# And finally JOIN the two together.
SL_joined_df <- full_join(SL_App_df, SL_Sav_df)

# Remove NA and NaN values, because these will cause the model to break.
SL_joined_clean <- SL_joined_df %>%
  filter(!is.na(PAR_surface)) %>% # removes NAs
  filter(!is.nan(PAR_surface)) # removes NaNs

# And create a list form as well.
SL_joined_list <- split(SL_joined_clean, SL_joined_clean$site_name)

## Save both formats, df and list.

saveRDS(SL_joined_df, "data_working/StreamLight_daily_df_192sites.rds")
saveRDS(SL_joined_list, "data_working/StreamLight_daily_list_192sites.rds")

# End of script.
