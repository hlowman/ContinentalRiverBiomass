## Resilience of Stream Productivity to Disturbance
## August 31, 2022
## Heili Lowman

# The following set of scripts will walk through the steps necessary to
# prep and send data to Teton as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_raw" and "data_working" folders have been ignored
# using git.ignore, so links to the raw data sets are provided in the step1
# file. If you are accessing the code via GitHub, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running
# the code.

# NOTE: This script is specifically designed to compared the Ricker model's
# calculations of "c" to critical discharge calculated based on
# geomorphic data provided by Jud Harvey.

#### Setup ####

## Load packages
lapply(c("tidyverse", "lubridate", "here", 
         "data.table", "pipeR", "patchwork"), require, character.only=T)

#### Data Import ####

# Load site-level info (hypoxia and Appling datasets).
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
site <- fread("data_raw/site_data.tsv")

# Load site-level parameters dataset with all iterations included.
dat_out <- readRDS("data_working/teton_190rivers_model_site_params_all_iterations_082422.rds")

# Load in original data fed into the model.
dat_in <- readRDS("data_working/df_190sites_10yrQnorm_allSL.rds")

#### Data Selection ####

# Geomorphic data is only available for the following sites:
# nwis_05451210: South Fork Iowa River, New Providence, IA
# nwis_05524500: Iroquois River, Foresman, IN
# nwis_05579630: Kickapoo Creek, Bloomington, IL
# nwis_01645704: Difficult Run, Fairfax, VA

# Filtering by those sites
my_list <- c("nwis_05451210","nwis_05524500",
             "nwis_05579630","nwis_01645704")

dat_out_4 <- dat_out %>%
  filter(site_name %in% my_list)

# And calculate normalized c statistic values
c_summary <- dat_out_4 %>%
  group_by(site_name) %>%
  summarize(mean_c = mean(c),
            med_c = median(c),
            q10_c = quantile(c, c(.10)),
            q90_c = quantile(c, c(.90))) %>%
  ungroup()

# Identify ten year flood values to back calculate c in cfs.
Q10_4 <- dat_in %>%
  filter(site_name %in% my_list) %>%
  distinct(site_name,RI_10yr_Q,RI_10yr_Q_cms)

# Join these datasets.
c_Q10_4 <- left_join(c_summary, Q10_4) # joined by "site_name"

# Calculate "c" in cubic meters per second.
c_Q10_4 <- c_Q10_4 %>%
  mutate(mean_c_cms = mean_c * RI_10yr_Q_cms,
         q50_c_cms = med_c * RI_10yr_Q_cms,
         q10_c_cms = q10_c * RI_10yr_Q_cms,
         q90_c_cms = q90_c * RI_10yr_Q_cms)

# Convert to cfs.
c_Q10_4 <- c_Q10_4 %>%
  mutate(mean_c_cfs = mean_c_cms * 35.3147,
         q50_c_cfs = q50_c_cms * 35.3147,
         q10_c_cfs = q10_c_cms * 35.3147,
         q90_c_cfs = q90_c_cms * 35.3147)

# Trim and export data.
c_Q10_4_trim <- c_Q10_4 %>%
  select(site_name, mean_c_cfs, q10_c_cfs, q50_c_cfs, q90_c_cfs)

write_csv(c_Q10_4_trim, 
          "data_working/critical_discharge_cfs_4_sites_090122.csv")

# End of script.
