# Data QAQC
# April 28, 2021
# Heili Lowman

# This script will be use for preliminary filtering
# of the Appling et al. 2018 dataset (citation below)
# to create a dataset of usable metabolism estimates
# based on the following criteria:

# (1) Sites should have a minimum of 275 days of data/year.

# (2) Sites should have data gaps no longer than 14 days.

# (3) Rhat values for GPP estimates should be less than 1.05.

# Some code has been adapted from Joanna Blaszczak's NWIS_RiverSelection
# script located at: 
# https://github.com/jrblaszczak/RiverBiomass/blob/main/code/NWIS_RiverSelection.R


# Setup -------------------------------------------------------------------

# Load packages.
library(tidyverse)
library(lubridate)
library(here)

# Double check working directory
here()

# Load datasets.
# Site information from A. Appling
sitedat <- read_tsv(here("data_raw", 
                         "site_data.tsv"))
# Model predictions from A. Appling
dailydat <- read_tsv(here("data_raw", 
                          "daily_predictions.tsv"))
# Site information from L. Koenig
hypoxdat <- read_csv(here("data_raw",
                          "GRDO_GEE_HA_NHD_2021_02_07.csv")) %>%
  filter(DB_Source == "PC") # filter only for Appling data


# Tidy --------------------------------------------------------------------

# Join site info datasets by siteID
sitesjoin <- left_join(sitedat, hypoxdat, by = c("site_name" = "SiteID"))
# Checked to make sure all 602 original sites were preserved.

# Export newly joined dataset to working file directory
saveRDS(sitesjoin, file = "data_working/sitesjoin.rds")


# Filter by Date ----------------------------------------------------------

# First, filter for amount of data using the metabolism dataset

# Format dates
dailydat <- dailydat %>%
  mutate(doy = yday(date), # create day of year column
         year = year(date)) # create year column

# Identify how many days for which there is data
# Count days per year
dat_per_year <- dailydat %>%
  group_by(site_name, year) %>%
  count() %>%
  ungroup()

# Create new dataset filtering by minimum of 275 days
minday275 <- dat_per_year %>%
  filter(n > 275)

# Identify max day gap per year
# Create new dataset with column added for daily gaps
gap_per_year <- dailydat %>%
  group_by(site_name, year) %>%
  mutate(gap = doy - lag(doy, default = doy[1])) %>% 
  # lag() for finding previous values
  ungroup()

# Create new dataset filtering by maximum daily gap of 14 days
maxgap14 <- gap_per_year %>%
  group_by(site_name, year) %>% # group by year
  summarize_at(.vars = "gap", .funs = max) %>% 
  ungroup() %>%
  filter(gap <= 14)

# Merge with days per year dataset
gaps_and_days <- merge(maxgap14, minday275, by = c("site_name", "year"))

# Examine how many unique sites this is
length(levels(factor(gaps_and_days$site_name))) # 203

# Filter by Data ----------------------------------------------------------

# Then, filter for quality of data using the metabolism dataset

minRhat105 <- dailydat %>%
  filter(GPP.Rhat < 1.05)

# Examine how many unique sites this is
length(levels(factor(minRhat105$site_name))) # 355

# Export Filtered Data ----------------------------------------------------

# Join filtered datasets
filterdat <- inner_join(gaps_and_days, minRhat105, by = c("site_name", "year"))

# Examine how many unique sites this is
length(levels(factor(filterdat$site_name))) # 202

# Export newly joined dataset to working file directory
saveRDS(filterdat, file = "data_working/filterdat.rds")

# End of script.
