## Resilience of Stream Productivity to Disturbance
## July 25, 2022
## Heili Lowman

# The following set of scripts will walk through the steps necessary to
# prep and send data to Teton as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_raw" and "data_working" folders have been ignored
# using git.ignore, so links to the raw data sets are provided below.
# If you are accessing the code via GitHub, these will need to be downloaded 
# and added to a folder of the appropriate name prior to running the code.

############################
## Setup
############################

lapply(c("tidyverse", "plyr", "cowplot","lubridate",
         "data.table","patchwork", "here"), require, character.only=T)

# Load in datasets created in NWIS_RiverSelection script.
NWIS_207sites_subset <- readRDS("data_working/NWIS_207sites_subset.rds")
NWIS_207sitesinfo_subset <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")
NWIS_207sites_Ndays <- readRDS("data_working/NWIS_207sites_Ndays.rds")

# Load in dataset created in StreamLight_data_assembly script.
SL_85sites <- readRDS("data_working/StreamLight_daily_85site_list.rds")

## Sites of interest that have StreamLight data
sites <- NWIS_207sitesinfo_subset$site_name
sites_files <- rep(NA, length(sites))
for(i in 1:length(sites)){
  sites_files[i] <- paste(sites[i],"_input_output.txt", sep="")
}

# need to identify the intersection of files with available light data
filenames <- intersect(sites_files, list.files(here("data_raw/individual_files")))

## Import StreamLight
SL <- ldply(filenames, function(filename) {
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
## Subset and evaluate NA
###########################
SF_df <- ldply(SL_daily, data.frame)

## add site_name
SF_df$site_name <- substr(SF_df$.id, 1, nchar(SF_df$.id)-17)

site_subset <- SF_df

site_subset_split <- split(site_subset, site_subset$.id)

# Check for NA values.
lapply(site_subset_split, function(x) sum(is.na(x$PAR_surface)))
lapply(site_subset_split, function(x) sum(is.na(x$PAR_turb)))

## Save

saveRDS(site_subset, "data_working/StreamLight_daily_207sites.rds")

# End of script.
