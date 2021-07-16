## 34 rivers data source
## July 13, 2021
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to practice fitting the Ricker model to ~30 site-years of data.
# The compiled data will then be run by sending the job to Teton (UWyoming).

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

# changed all directories to match folder in modelscape directory on Teton
data <- readRDS("NWIS_34sites_subset.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("NWIS_34sitesinfo_subset.rds")

# Read in StreamLight processed data (Savoy)
SL <- readRDS("StreamLight_daily_34sites.rds")
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

## visualize
# ggplot(data, aes(date, GPP)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(~site_name,scales = "free_x")

## split list by ID
l <- split(data, data$site_name)

## For sites without StreamLight data, set PAR_surface to light
# first create a function to map across the list
swap_light <- function(df){
  df %>%
    mutate(PAR_new = case_when(is.na(PAR_surface) == TRUE ~ light,
                               TRUE ~ PAR_surface))
}

# test the above function on a single dataframe that doesn't have light
# test_df <- l$nwis_01656903
# test_output <- swap_light(test_df) # Yay! :)
# and for one that does have light...
# test2_df <- l$nwis_01608500
# test2_output <- swap_light(test2_df) # Double Yay! :)

l <- map(l, swap_light)

rel_LQT <- function(x){
  x$light_rel <- x$PAR_new/max(x$PAR_new)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)

  x<-x[order(x$date),]
  return(x)
}

dat <- lapply(l, function(x) rel_LQT(x))
df <- dat

rm(data, l, SL, data_siteyears, dat)

# End of script.
