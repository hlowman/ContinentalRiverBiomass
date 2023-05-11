## Recovery of Stream Productivity following Disturbance
## Originally created: December 1, 2022
## Heili Lowman

#### READ ME ####

# The following set of scripts will walk through the steps necessary to
# prep and send data to Beartooth as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_XX" folders have been ignored using git.ignore.
# Links to the raw data sets are provided in the 01_NWIS_RiverSelection.R file.

# If you are accessing the code via GitHub, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running the 
# code below.

# The following script will calculate the number of times each site
# exceeds what we have estimated as the critical disturbance threshold (Qc).

#### Setup ####

# Load packages.
lapply(c("lubridate","tidyverse", "here", "viridis",
         "reshape2","ggExtra","patchwork"), require, character.only=T)

# Load datasets.

# First, the raw data loaded into the model.
dat_in <- readRDS("data_working/list_181sites_Qmaxnorm_SavoySL.rds")

# And then a dataset containing Qc values that have been converted from
# c values estimated by our model (all, not filtered by diagnostics).
dat_Qc <- readRDS("data_working/QcQ2_unfiltered_181sites_051023.rds")

# Also also, hypoxia dataset for additional info re: slope
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")

#### Test Workflow ####

# Now, to create a practice workflow to later turn into a function to
# quantify exceedances of Qc. Going to first pick a site where it actually
# happens.

test_df <- dat_in$nwis_01400500

# First visualize the data to see how many times it happens
ggplot(test_df, aes(x = Date, y = Q)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 67.4526324) +
  labs(x = "Date",
       y = "Discharge") +
  theme_bw()
# Ok, so roughly 20 times.

# For ease, making the original list into a df.
dat_in_df <- map_df(dat_in, ~as.data.frame(.x), .id="site_name")

# Trim down to columns of interest for Qc dataset
dat_Qc_trim <- dat_Qc %>%
  dplyr::select(site_name, Qc)

# Need to add Qc values for each site to the input data
dat_in_Qc <- inner_join(dat_in_df, dat_Qc_trim)

# And make back into a list to apply the function below.
## split list by ID
dat_in_Qc_l <- split(dat_in_Qc, dat_in_Qc$site_name)

#### Exceedance Workflow ####

# Modifying my event delineation function from step3 of this workflow.

# Event delineation function:
# loop over the separate time sequences for a given site
exceedFUN <- function(d){
  
  # calculate the difference from one day to the next - YES/NO
    
    d <- d %>%
      mutate(exceedance = case_when(Q > Qc ~ "YES", TRUE ~ "NO"))
  
  # delineate sequenced time frames based on day to day differences
    
    d <- d %>% # mark all yesses preceded by nos
      mutate(sequence = case_when(exceedance == "YES" & 
                               lag(exceedance) == "NO" ~ 1,
                             TRUE ~ 0)) %>%
      # and count total days of exceedance
      mutate(total = case_when(exceedance == "YES" ~ 1,
                               TRUE ~ 0))
  
  return(d)
  
}

# apply event delineation function
dat_in_exc <- lapply(dat_in_Qc_l, function(x) exceedFUN(x))

# spot checking a site to see how it looks
site1 <- dat_in_exc$nwis_02203900

# Make the list back into a df
dat_in_exc_df <- map_df(dat_in_exc, ~as.data.frame(.x), .id="site_name")

# And summarize total exceedance events and exceedance days as well as
# total days in a record.
dat_exceed_sum <- dat_in_exc_df %>%
  group_by(site_name) %>%
  summarize(total_exc_events = sum(sequence),
            total_exc_days = sum(total),
            total_days = n()) %>%
  ungroup() %>%
  # calculate # of exceedance days relative to length of dataset
  mutate(total_exc_rel_d = total_exc_days/total_days)

# Standardize by years on record for each site.
dat_years <- dat_in_df %>%
  group_by(site_name) %>%
  summarize(years = n_distinct(year)) %>%
  ungroup()

# Join with accrual dataset.
dat_exc_y <- inner_join(dat_exceed_sum, dat_years)

# and make a column with # of exceedances normalized by # of years on record.
dat_exc_y <- dat_exc_y %>%
  mutate(exc_y = total_exc_events/years)

# And export the exceedance data for use in the linear models.
saveRDS(dat_exc_y, "data_working/Qc_exceedances_181sites_051123.rds")

# End of script.
