## Recovery of Stream Productivity following Disturbance
## Originally created: October 12, 2022
## Heili Lowman

#### READ ME ####

# The following set of scripts will walk through the steps necessary to
# prep and send data to the University of Wyoming's Beartooth High Performance
# Computing Cluster as well as process the model outputs.

# Much of this code has been adapted from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_XX" folders have been ignored using git.ignore.
# Links to the raw data sets are provided in the 01_NWIS_RiverSelection.R 
# and 02_StreamLight_processing files.

# If you are accessing the code via GitHub or Zenodo, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running the 
# code below.

# This script compiles GPP, light, and discharge data to use when fitting
# the latent biomass models.

#### Setup ####

## Load packages
lapply(c("tidyverse", "cowplot",
         "lubridate", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", "here"), require, character.only=T)

#### Data Import & Processing #####

# Import stream daily GPP data
data <- readRDS("data_working/NWIS_198sites_050423.rds")
data$Date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

# Import stream metadata info
site_info <- readRDS("data_working/NWIS_198sitesinfo_050423.rds")

# Read in StreamLight data (Savoy)
SL <- readRDS("data_working/StreamLight_daily_df_182sites_050423.rds")

## Join GPP data and StreamLight
data_join <- left_join(data, SL, by=c("site_name", "Date"))
# double-check, 213,519 observations still

# But need to again remove days for which we don't have light estimates.
data_join_tidylight <- data_join %>%
  filter(!is.na(PAR_surface))
# Now down to 194,100 observations

## How many days of data per site per year
data_join_tidylight$year <- year(data_join_tidylight$date)

data_siteyears <- data_join_tidylight %>%
  dplyr::count(site_year) # 684 unique site-years

sum(data_siteyears$n) # 194,100 days ~ 532 years

## Set any GPP < 0 to a small value between 0.05 to 0.13 g O2 m-2 d-1
# Doing this a second time, because the first time this was done in step 1
# it was only to calculate mean/max GPP and it was saved as a separate
# column, "GPP_temp".
set.seed(148)
data_join_tidylight[which(data_join_tidylight$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
# Assuming a normal distribution: mean - 1.96*SD = lower CI (2.5%tile)
data_join_tidylight$GPP_sd <- (((data_join_tidylight$GPP.upper - data_join_tidylight$GPP)/1.96) + 
                       ((data_join_tidylight$GPP.lower - data_join_tidylight$GPP)/-1.96))/2

# Addition of time-series delineation column using seqFUN function created below
# Borrowing this code from code/teton_moresites/Data_Availability_figures.R

# Made sure the light data that needed to be removed (due to NAs and NaNs) was
# taken care of above, since this would affect the gaps in data.

# Shortening the name for ease of use.
data_tidy <- data_join_tidylight

# create placeholder columns for day-to-day differences and sequences
data_tidy <- data_tidy %>%
  mutate(diff_time = 0, seq = 1, new_e = 0)

# Addition of site numbers for model for loops.
sites <- unique(data_tidy$site_name)
index <- seq(1:181) 
sites_index <- as.data.frame(cbind(index, sites))

data_tidy <- left_join(data_tidy, sites_index, by = c("site_name" = "sites"))

# Addition of day numbers for further model indexing.
data_tidy <- left_join(data_tidy, data_siteyears) %>%
  dplyr::rename(Ndays = n)
# Looks like no sites have dropped below 30 days worth of records following
# removal of NAs for light.

# Due to date issues below, forcing the dataframe to present dates in ascending order
data_tidy <- data_tidy %>%
  group_by(site_name) %>%
  arrange(date) %>% # ascending is the default
  ungroup()

## Re-identifying the max day gap PER YEAR within a given site
# to be sure no remaining site_years need to be filtered out following
# light data removal.
gap_per_year <- data_tidy %>%
  group_by(site_name, year) %>%
  dplyr::mutate(gap = doy - lag(doy, default=doy[1], order_by = doy)) %>%
  ungroup()
# The code above really struggles between the dplyr/plyr packages it seems.
# So, ALWAYS double check the output above to be sure there aren't any gaps
# that are negative, as that indicates the data isn't being grouped properly.

maxgap <- gap_per_year %>%
  group_by(site_name, year) %>%
  summarize_at(.vars = "gap", .funs = max) %>%
  ungroup()
# Ok - we're in the clear, all gaps are still 14 days or less.

## split list by ID
l <- split(data_tidy, data_tidy$site_name)

# Event delineation function:
# loop over the separate time sequences for a given site
seqFUN <- function(d){
  
  # calculate the difference from one day to the next
  for(i in 2:nrow(d)){
    d$diff_time[i] = difftime(time1 = d$date[i], time2 = d$date[(i-1)],
                              units = "days")
  }
  
  # delineate sequenced time frames based on day to day differences
  # anything less than 14 day gaps is permissible
  for(i in 2:nrow(d)){
    if(d$diff_time[i] < 14){
      d$seq[i] = d$seq[(i-1)]
    } else {
      d$seq[i] = d$seq[(i-1)] + 1
    }
  }
  
  # add column to delineate changes in events
  for(i in 2:nrow(d)){
    d$new_e[i] = d$seq[i] - d$seq[i-1]
  }
  
  return(d)
  
}

# apply event delineation function
l <- lapply(l, function(x) seqFUN(x))

# Do a spot check on a few sites to be sure (a) delineation function works and 
# (b) no NAs or NaNs in light values. - Looks good!
View(l$nwis_02208493)
View(l$nwis_02203700)
View(l$nwis_04142000)

# and a quick histogram of all light values
hist(data_tidy$PAR_surface)

# Create function for relative light as well as
# standardized discharge columns relative to !Q max! during a given record
rel_LQ <- function(x){
  x$light_rel <- x$PAR_surface/max(x$PAR_surface, na.rm = TRUE)

  x$Q_rel <- x$Q/max(x$Q, na.rm = TRUE)

  x<-x[order(x$date),]
  return(x)
}

# apply relativizing/standardizing function
dat_spring23 <- lapply(l, function(x) rel_LQ(x))

# Turn back into a dataframe for final filtering steps.
dat_spring23_df <- plyr::ldply(dat_spring23, data.frame)

# Double check to be sure how many sites remain.
length(unique(dat_spring23_df$site_name)) # 181

# Exporting datasets, both formats
saveRDS(dat_spring23, "data_working/list_181sites_Qmaxnorm_SavoySL.rds")
saveRDS(dat_spring23_df, "data_working/df_181sites_Qmaxnorm_SavoySL.rds")

# Also breaking the list up into 4 parts for parallel jobs.
# this makes model fitting faster, troubleshooting easier, and 
# load-in of the data outputs quicker later.
dat_list_pt1 <- dat_spring23[1:45]
dat_list_pt2 <- dat_spring23[46:90]
dat_list_pt3 <- dat_spring23[91:135]
dat_list_pt4 <- dat_spring23[136:181]

# Exporting these four lists
saveRDS(dat_list_pt1, "data_working/list_181sites_Qmaxnorm_SavoySL_pt1.rds")
saveRDS(dat_list_pt2, "data_working/list_181sites_Qmaxnorm_SavoySL_pt2.rds")
saveRDS(dat_list_pt3, "data_working/list_181sites_Qmaxnorm_SavoySL_pt3.rds")
saveRDS(dat_list_pt4, "data_working/list_181sites_Qmaxnorm_SavoySL_pt4.rds")

# End of script.
