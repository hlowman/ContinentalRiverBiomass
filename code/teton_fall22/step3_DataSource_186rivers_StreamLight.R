## Resilience of Stream Productivity to Disturbance
## October 12, 2022
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

#### Setup ####

## Load packages
lapply(c("tidyverse", "cowplot",
         "lubridate", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", "here"), require, character.only=T)

#### Data Import & Processing #####

# Import stream GPP data
data <- readRDS("data_working/NWIS_198sites_101222.rds")
data$Date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

# Import stream info
site_info <- readRDS("data_working/NWIS_198sitesinfo_101222.rds")

# Read in StreamLight processed data (Appling+Savoy)
SL <- readRDS("data_working/StreamLight_daily_df_186sites.rds")

## Join GPP data and StreamLight
data_join <- left_join(data, SL, by=c("site_name", "Date"))

# But need to again remove days for which we don't have light estimates
# otherwise the model freaks out.
data_join_tidylight <- data_join %>%
  filter(!is.na(PAR_surface))

## How many days of data per site per year
data_join_tidylight$year <- year(data_join_tidylight$date)

data_siteyears <- data_join_tidylight %>%
  count(site_name, year) # 691 unique site-years

sum(data_siteyears$n) # 196,963 days ~ 540 years

## Set any GPP < 0 to a small value between 0.05 to 0.13 g O2 m-2 d-1
# Doing this a second time, because the first time this was done in step 1
# it was only to calculate mean/max GPP and it was saved as a separate
# column, "GPP_temp".
data_join_tidylight[which(data_join_tidylight$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
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
index <- seq(1:182) # Be sure to change this when having changed filters!
sites_index <- as.data.frame(cbind(index, sites))

data_tidy <- left_join(data_tidy, sites_index, by = c("site_name" = "sites"))

# Addition of day numbers for further model indexing.
data_tidy <- left_join(data_tidy, data_siteyears) %>%
  rename(Ndays = n)

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

# Create function for relative light and temperature as well as
# standardized discharge columns relative to !Q max! during a given record
rel_LTQ <- function(x){
  x$light_rel <- x$PAR_surface/max(x$PAR_surface, na.rm = TRUE)
  x$temp_rel <- x$temp/max(x$temp, na.rm = TRUE)
  x$Q_rel <- x$Q/max(x$Q, na.rm = TRUE)

  x<-x[order(x$date),]
  return(x)
}

# apply relativizing/standardizing function
dat_fall22 <- lapply(l, function(x) rel_LTQ(x))

# Turn back into a dataframe for final filtering steps.
dat_fall22_df <- plyr::ldply(dat_fall22, data.frame)

# Double check to be sure how many sites remain.
length(unique(dat_fall22_df$site_name)) #182

# Exporting datasets, both formats
saveRDS(dat_fall22, "data_working/list_182sites_Qmaxnorm_allSL.rds")
saveRDS(dat_fall22_df, "data_working/df_182sites_Qmaxnorm_allSL.rds")

# Also breaking the list up into 3 parts for parallel jobs.
dat_list_pt1 <- dat_fall22[1:60]
dat_list_pt2 <- dat_fall22[61:120]
dat_list_pt3 <- dat_fall22[121:182]

# Exporting these three lists
saveRDS(dat_list_pt1, "data_working/list_182sites_Qmaxnorm_allSL_pt1.rds")
saveRDS(dat_list_pt2, "data_working/list_182sites_Qmaxnorm_allSL_pt2.rds")
saveRDS(dat_list_pt3, "data_working/list_182sites_Qmaxnorm_allSL_pt3.rds")

# End of script.
