## 34 rivers data source
## Step THREE in Metabolism Modeling Workflow
## November 16, 2021
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to practice fitting the Ricker model to a series of selected sites.

# Specifically, the workflow in this folder will work on pooling data,
# working across gaps in time series as well as models w/ and w/o the P term.

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
data <- readRDS("data_working/NWIS_pooling_3sites_subset.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("data_working/NWIS_pooling_3sitesinfo_subset.rds")

## Removed code adding Stream Light data because there was none
## for these three sites.

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

## Add new column called PAR_new to keep naming convention consistent
## despite lack of StreamLight data
data <- data %>%
  mutate(PAR_new = light)

## split list by ID
l <- split(data, data$site_name)

## Removed function swapping in StreamLight data since there is none.

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

# Exporting dataset
saveRDS(df, "data_working/df_pooling_3sites.rds") # data for pooling test runs

# End of script.
