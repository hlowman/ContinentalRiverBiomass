## 3 rivers data source
## Created: October 13, 2021
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to practice fitting the Ricker model to 3 sites of data.

# I've commented out those steps that I feel I don't
# need to perform. I've also changed the filepaths to match my
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
data <- readRDS("data_working/NWIS_3sites_subset.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("data_working/NWIS_3sitesinfo_subset.rds")

# Read in StreamLight processed data (Savoy)
SL <- readRDS("data_working/StreamLight_daily_3sites.rds")
colnames(SL)[colnames(SL) == "Date"] <- "date"

## Join data and StreamLight
data <- left_join(data, SL, by=c("site_name", "date"))

## How many days of data per site per year
data$year <- year(data$date)

data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()

## Create additional gappy dataset
y1 <- c(2008, 2009, 2011, 2012, 2013, 2016) # removing 2010, 2014-15
y2 <- c(2008, 2009, 2010, 2013, 2015, 2016) # removing 2012, 2014
y3 <- c(2008, 2010, 2012, 2013, 2014, 2015) # removing 2016

data_w_gaps <- data %>%
  filter(site_name == "nwis_02266300" & year %in% y1 |
           site_name == "nwis_14206950" & year %in% y2 |
           site_name == "nwis_01608500" & year %in% y3)

## Set any GPP < 0 to a small value between 0.05 to 0.13 g O2 m-2 d-1
data[which(data$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)
data_w_gaps[which(data_w_gaps$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2
data_w_gaps$GPP_sd <- (((data_w_gaps$GPP.upper - data_w_gaps$GPP)/1.96) + ((data_w_gaps$GPP.lower - data_w_gaps$GPP)/-1.96))/2

## split list by ID
l <- split(data, data$site_name)
l_gap <- split(data_w_gaps, data_w_gaps$site_name)
## and for temporal comparison, year-ID of site that performs well
one_site <- data %>% filter(site_name == "nwis_01608500")
l_yrs <- split(one_site, one_site$site_year)

## For sites without StreamLight data, set PAR_surface to light
# first create a function to map across the list
swap_light <- function(df){
  df %>%
    mutate(PAR_new = case_when(is.na(PAR_surface) == TRUE ~ light,
                               TRUE ~ PAR_surface))
}

l <- map(l, swap_light)
l_gap <- map(l_gap, swap_light)
l_yrs <- map(l_yrs, swap_light)

rel_LQT <- function(x){
  x$light_rel <- x$PAR_new/max(x$PAR_new)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)

  x<-x[order(x$date),]
  return(x)
}

dat <- lapply(l, function(x) rel_LQT(x))
df <- dat

dat_gap <- lapply(l_gap, function(x) rel_LQT(x))
df_gap <- dat_gap

dat_yrs <- lapply(l_yrs, function(x) rel_LQT(x))
df_yrs <- dat_yrs

#rm(data, l, SL, data_siteyears, dat)
#rm(data_gap, l_gap, SL, data_siteyears, dat_gap)

# Exporting dataset
saveRDS(df, "data_working/df_3sites.rds") # first dataset for test run
saveRDS(df_gap, "data_working/df_3sites_gappy.rds") # second dataset for test run
saveRDS(df_yrs, "data_working/df_1site_allyrs.rds") # third dataset for test run - separated by year rather than site name

# End of script.
