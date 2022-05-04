## 207 rivers data source
## Step THREE in Metabolism Modeling Workflow
## April 20, 2022
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to fit the Ricker model to 206 sites for my JASM poster presentation.

# I've commented out those steps that I feel, for the time being, I don't
# need to perform, and I've changed the appropriate filepaths to match my
# repository structure.

## Load packages
lapply(c("tidyverse", "cowplot",
         "lubridate", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", "here"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################

# changed all directories back to match my folders
data <- readRDS("data_working/NWIS_207sites_subset.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")

# Read in StreamLight processed data (Savoy)
SL <- readRDS("data_working/StreamLight_daily_207sites.rds")
colnames(SL)[colnames(SL) == "Date"] <- "date"

# Read in 10 year flood data
q10 <- read_csv("data_working/RI_10yr_flood_206riv.csv")

## Join data and StreamLight
data <- left_join(data, SL, by=c("site_name", "date"))

## Join data and q10
data <- left_join(data, q10, by = c("site_name"))

## How many days of data per site per year
data$year <- year(data$date)

data_siteyears <- data %>%
  count(site_name, year)

sum(data_siteyears$n) # 219,871 days = 602.39 years

## DONT FORGET THIS:
## Set any GPP < 0 to a small value between 0.05 to 0.13 g O2 m-2 d-1
data[which(data$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2

# Addition of time-series delineation column using seqFUN function created below
# Borrowing this code from code/teton_moresites/Data_Availability_figures.R
# create placeholder columns for day to day differences and sequences
data <- data %>%
  mutate(diff_time = 0, seq = 1, new_e = 0)

# Addition of site numbers for model for loops.
sites <- unique(data$site_name)
index <- seq(1:207)
sites_index <- as.data.frame(cbind(index, sites))

data <- left_join(data, sites_index, by = c("site_name" = "sites"))

# Addition of day numbers for further model indexing.
data <- left_join(data, data_siteyears) %>%
  rename(Ndays = n)

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
l <- lapply(l, function(x) seqFUN(x))

## For sites without StreamLight data, set PAR_surface to light
# first create a function to map across the list
swap_light <- function(df){
  df %>%
    mutate(PAR_new = case_when(is.na(PAR_surface) == TRUE ~ light,
                               TRUE ~ PAR_surface))
}

l <- map(l, swap_light)

# Create function for relative light and temperature as well as
# standardized discharge columns
rel_LQT <- function(x){
  x$light_rel <- x$PAR_new/max(x$PAR_new)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/x$RI_10yr_Q_cms # standardizing by 10 year flood in cubic meters/second

  x<-x[order(x$date),]
  return(x)
}

# apply relativizing/standardizing function
dat <- lapply(l, function(x) rel_LQT(x))

#Note - no discharge data for nwis 03293500 so removing it
# trim and rename dataset
df <- dat[-89]

# Exporting dataset
saveRDS(df, "data_working/df_206sites_10yrQnorm.rds") # data for random effects testing runs

# End of script.
