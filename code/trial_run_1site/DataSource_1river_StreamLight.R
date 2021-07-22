## 1 river data source
## Created: June 23, 2021
## Heili Lowman

# I'll be modifying some of Joanna's code from the RiverBiomass repository
# to practice fitting the Ricker model to 1 year of data at a "good" site.

# I've commented out those steps that I feel, for the time being, I don't
# need to perform, and I've changed the appropriate filepaths to match the
# file structure on Teton.

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", "here"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################
data <- readRDS("/project/modelscape/users/hlowman/jobscripts/1site_test/NWIS_1site_subset_good_2015.rds")

data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("/project/modelscape/users/hlowman/jobscripts/1site_test/NWIS_1siteinfo_subset_good_2015.rds")

# Read in StreamLight processed data (Savoy)
SL <- readRDS("/project/modelscape/users/hlowman/jobscripts/1site_test/StreamLight_daily_1riv_good_2015.rds")

colnames(SL)[colnames(SL) == "Date"] <- "date"

## Join data and StreamLight
data <- left_join(data, SL, by=c("site_name", "date"))

## Change river names to short names
site_info$short_name <- revalue(as.character(site_info$site_name), replace = c("nwis_01608500"="South Branch Potomac River, WV"))

## How many days of data per site per year
data$year <- year(data$date)

data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()

## Select the first of two years
data <- rbind(data[which(data$site_name == "nwis_01608500" & data$year %in% c(2015)),])

## Set any GPP < 0 to a small value between 0.05 to 0.13 g O2 m-2 d-1
data[which(data$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2

## visualize
# ggplot(data, aes(date, GPP))+
#   geom_point()+geom_line()+
#   facet_wrap(~site_name,scales = "free_x")

## split list by ID
l <- split(data, data$site_name)

rel_LQT <- function(x){
  x$light_rel <- x$PAR_surface/max(x$PAR_surface)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)

  x<-x[order(x$date),]
  return(x)
}

dat <- lapply(l, function(x) rel_LQT(x))
df <- dat

rm(data, l, SL, data_siteyears, dat)

# End of script.
