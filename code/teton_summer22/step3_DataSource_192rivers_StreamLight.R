## Resilience of Stream Productivity to Disturbance
## July 29, 2022
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

## Load packages
lapply(c("tidyverse", "cowplot",
         "lubridate", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", "here"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################

# Import stream GPP data
data <- readRDS("data_working/NWIS_207sites_subset2022.rds")
data$Date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

# Import stream info
site_info <- readRDS("data_working/NWIS_207sitesinfo_subset2022.rds")

# Read in StreamLight processed data (Appling+Savoy)
SL <- readRDS("data_working/StreamLight_daily_df_192sites.rds")

# Read in 10 year flood data
q10 <- read_csv("data_working/RI_10yr_flood_206riv.csv")

## Join GPP data and StreamLight
data_join <- left_join(data, SL, by=c("site_name", "Date"))

## Join this dataset and q10
data_join <- left_join(data_join, q10, by = c("site_name"))

## How many days of data per site per year
data_join$year <- year(data_join$date)

data_siteyears <- data_join %>%
  count(site_name, year) # 784 unique site-years

sum(data_siteyears$n) # 219,871 days = 602.39 years

## Set any GPP < 0 to a small value between 0.05 to 0.13 g O2 m-2 d-1
data_join[which(data_join$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data_join$GPP_sd <- (((data_join$GPP.upper - data_join$GPP)/1.96) + 
                       ((data_join$GPP.lower - data_join$GPP)/-1.96))/2

# Addition of time-series delineation column using seqFUN function created below
# Borrowing this code from code/teton_moresites/Data_Availability_figures.R

# create placeholder columns for day to day differences and sequences
data_join <- data_join %>%
  mutate(diff_time = 0, seq = 1, new_e = 0)

# Addition of site numbers for model for loops.
sites <- unique(data_join$site_name)
index <- seq(1:207)
sites_index <- as.data.frame(cbind(index, sites))

data_join <- left_join(data_join, sites_index, by = c("site_name" = "sites"))

# Addition of day numbers for further model indexing.
data_join <- left_join(data_join, data_siteyears) %>%
  rename(Ndays = n)

## split list by ID
l <- split(data_join, data_join$site_name)

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

# Create function for relative light and temperature as well as
# standardized discharge columns relative to 10yr flood
rel_LTQ <- function(x){
  x$light_rel <- x$PAR_surface/max(x$PAR_surface, na.rm = TRUE)
  x$temp_rel <- x$temp/max(x$temp, na.rm = TRUE)
  x$Q_rel <- x$Q/x$RI_10yr_Q_cms # standardizing by 10 year flood in cubic meters/second

  x<-x[order(x$date),]
  return(x)
}

# apply relativizing/standardizing function
dat22 <- lapply(l, function(x) rel_LTQ(x))

# Turn back into a dataframe for final filtering steps.
dat22_df <- plyr::ldply(dat22, data.frame)

# No discharge data for nwis 03293500, so need to remove it.
# But can't remove just NA records, because this will capture other
# sites as well, so first need to take the mean RI_10yr_Q_cms at each site.
means <- dat22_df %>%
  group_by(site_name) %>%
  summarize(meanQ = mean(RI_10yr_Q_cms, na.rm = TRUE)) %>%
  ungroup()

means_na <- means %>%
  filter(is.na(meanQ))

# Filter by missing mean.
dat22_df1 <- dat22_df %>%
  filter(site_name != "nwis_03293500")

# Double check to be sure it was removed.
length(unique(dat22_df$site_name)) # 207
length(unique(dat22_df1$site_name)) # 206

# For sites without StreamLight data, remove them.
meansSL <- dat22_df1 %>%
  group_by(site_name) %>%
  summarize(meanSL = mean(PAR_surface, na.rm = TRUE)) %>%
  ungroup()

meansSL_na <- meansSL %>%
  filter(is.na(meanSL)) # 16 sites

# Make a list of these sites
to_remove <- meansSL_na$site_name

# And filter by this liste
dat22_df2 <- dat22_df1 %>%
  filter(!site_name %in% to_remove)

# Double check to be sure the sites that should have been removed
# were indeed removed.
length(unique(dat22_df1$site_name)) # 206
length(unique(dat22_df2$site_name)) # 190

# And did a few spot checks for sites that shouldn't be in there.

# Create the final dataset in list form as well.
dat22_list <- split(dat22_df2, dat22_df2$.id)

# Exporting datasets, both formats
saveRDS(dat22_df2, "data_working/df_190sites_10yrQnorm_allSL.rds")
saveRDS(dat22_list, "data_working/list_190sites_10yrQnorm_allSL.rds")

# End of script.
