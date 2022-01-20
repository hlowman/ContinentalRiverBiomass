## 4 rivers data source
## Step THREE in Metabolism Modeling Workflow
## January 20, 2022
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to practice fitting the Ricker model to a series of selected sites.

# Specifically, the workflow in this folder will work on 
#gaps in time series as well as models w/ and w/o the P term.

# I've commented out those steps that I feel, for the time being, I don't
# need to perform, and I've changed the appropriate filepaths to match my
# repository structure.

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", "here"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################

# changed all directories back to match my folders
data <- readRDS("data_working/NWIS_4sites_subset.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("data_working/NWIS_4sitesinfo_subset.rds")

# Read in StreamLight processed data (Savoy)
SL <- readRDS("data_working/StreamLight_daily_4sites.rds")
colnames(SL)[colnames(SL) == "Date"] <- "date"

## Join data and StreamLight
data <- left_join(data, SL, by=c("site_name", "date"))

## Removed code changing names to short names - see previous versions if 
## needed to add back in.

## How many days of data per site per year
data$year <- year(data$date)

data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()

## Using all data so removed code selecting only 2 years of data

## Set any GPP < 0 to a small value between 0.05 to 0.13 g O2 m-2 d-1
data[which(data$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2

## For sites without StreamLight data, set PAR_surface to light
# first create a function to map across the list
# swap_light <- function(df){
#   df %>%
#     mutate(PAR_new = case_when(is.na(PAR_surface) == TRUE ~ light,
#                                TRUE ~ PAR_surface))
# }
# 
# l <- map(l, swap_light)

# the above was giving me issues so i created the column
# before turning it into a list instead
data <- data %>%
  mutate(PAR_new = case_when(is.na(PAR_surface) == TRUE ~ light,
                             TRUE ~ PAR_surface))

# Addition of column for model selection EDIT THIS TO INCLUDE FILTER IN THE FUTURE
data <- data %>%
  mutate(p_remove = case_when(site_name %in% c("nwis_02266300", "nwis_03067510") ~ 1,
                              TRUE ~ 0))

# Addition of time-series delineation column using seqFUN function created below
# Borrowing this code from code/teton_moresites/Data_Availability_figures.R
# create placeholder columns for day to day differences and sequences
data <- data %>%
  mutate(diff_time = 0, seq = 1, new_e = 0)

## split list by ID
l <- split(data, data$site_name)

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
testing <- lapply(l, function(x) seqFUN(x))

# Create function for relative light and temperature as well as
# standardized discharge columns
rel_LQT <- function(x){
  x$light_rel <- x$PAR_new/max(x$PAR_new)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)

  x<-x[order(x$date),]
  return(x)
}

# apply relativizing/standardizing function
dat <- lapply(l, function(x) rel_LQT(x))

# rename dataset
df <- dat

#rm(data, l, SL, data_siteyears, dat)

# Exporting dataset
saveRDS(df, "data_working/df_4sites.rds") # data for pooling test runs

# End of script.
