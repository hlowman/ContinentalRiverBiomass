## Recovery of Stream Productivity following Disturbance
## Originally created: October 13, 2022
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

#### Setup ####

# Load necessary packages.
lapply(c("tidyverse", "lubridate", "data.table",
         "rstan","bayesplot","shinystan", "here",
         "GGally", "glmmTMB", "MuMIn", "effects",
         "DHARMa", "lme4", "multcomp", "patchwork",
         "calecopal", "viridis", "plotly", "ggbreak"), require, character.only=T)

# Load necessary datasets.
# Load site-level info (hypoxia and Appling datasets).
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
site <- fread("data_raw/site_data.tsv")

# Load in original list of data fed into the model.
dat_in <- readRDS("data_working/list_181sites_Qmaxnorm_SavoySL.rds")

# Load in model diagnostics for site-level parameters.
dat_diag <- readRDS("data_working/beartooth_181rivers_model_diags_050923.rds")

# Load in list containing all iterations of site-level parameters.
dat_out <- readRDS("data_working/beartooth_181rivers_model_params_all_iterations_050823.rds")

# Load 2yr flood data.
dat_2yrQ <- read_csv("data_working/RI_2yr_flood_180riv_050923.csv")

#### Data Prep ####

# Take list containing all input data and make into a df.
dat_in_df <- map_df(dat_in, ~as.data.frame(.x), .id="site_name")

# Take list containing all iterations of parameters and make into a df.
dat_out_df <- map_df(dat_out, ~as.data.frame(.x), .id="site_name")

#### Value filter for rmax ####

# Negative rmax values are not biologically reasonable, so I'm removing them. 
# Instead of calculating the mean, I'll be using median rmax values.
dat_out_rmed <- dat_out_df %>%
  group_by(site_name) %>%
  summarize(r_med = median(r)) %>%
  ungroup()

# And remove negative values.
dat_out_rmed_pos <- dat_out_rmed %>%
  filter(r_med > 0) # Removes 14 sites.

#### Rhat filter for rmax ####

# Before proceeding with the first step on my analyses, I will be filtering out 
# sites at which the model did not converge well for the rmax parameter.
# Sites with Rhat > 1.05 will not pass muster.

dat_diag_rfilter1 <- dat_diag %>%
  filter(parameter == "r") %>%
  filter(Rhat < 1.05) # 17 sites drop off

# Next, append the positive rmax values to the Rhat filter to remove
# appropriate sites.
dat_out_rmed_Rhat <- inner_join(dat_diag_rfilter1, dat_out_rmed_pos) 
# 152 sites remaining

# Export for future use.
saveRDS(dat_out_rmed_Rhat, "data_working/rmax_filtered_152sites_051023.rds")

#### Conversion of c to Qc ####

# Now, before performing any filtering for c values, I'm first going to calculate
# Qc:Q2yr so that I have these values for the manuscript table (since some don't
# pass c filtering steps and I don't want to have to back-track).

dat_out_cmed <- dat_out_df %>%
  group_by(site_name) %>%
  summarize(c_med = median(c)) %>%
  ungroup()

dat_diag_c <- dat_diag %>%
  filter(parameter == "c")

dat_out_cmed_diag <- left_join(dat_out_cmed, dat_diag_c)

# Now, convert normalized c values to typical discharge values.
dat_maxQ <- dat_in_df %>%
  group_by(site_name) %>%
  summarize(maxQ = max(Q, na.rm = TRUE)) %>%
  ungroup()

dat_out_cmed_maxQ <- left_join(dat_out_cmed_diag, dat_maxQ)

# convert both c and confidence interval values.
dat_out_cmed_maxQ$Qc <- dat_out_cmed_maxQ$c_med*dat_out_cmed_maxQ$maxQ
dat_out_cmed_maxQ$Qc2.5 <- dat_out_cmed_maxQ$`2.5%`*dat_out_cmed_maxQ$maxQ
dat_out_cmed_maxQ$Qc97.5 <- dat_out_cmed_maxQ$`97.5%`*dat_out_cmed_maxQ$maxQ

# Join with 2yr flood data.
dat_out_cmed_2yrQ <- left_join(dat_out_cmed_maxQ, dat_2yrQ)

# And determine Qc:Q2yr ratio values.
dat_out_cmed_2yrQ$Qc_Q2yr <- dat_out_cmed_2yrQ$Qc/dat_out_cmed_2yrQ$RI_2yr_Q_cms
dat_out_cmed_2yrQ$Qc_Q2yr2.5 <- dat_out_cmed_2yrQ$Qc2.5/dat_out_cmed_2yrQ$RI_2yr_Q_cms
dat_out_cmed_2yrQ$Qc_Q2yr97.5 <- dat_out_cmed_2yrQ$Qc97.5/dat_out_cmed_2yrQ$RI_2yr_Q_cms

# Export for future use.
saveRDS(dat_out_cmed_2yrQ, "data_working/QcQ2_unfiltered_181sites_051023.rds")

#### Value filter for c ####

# Negative c values are not biologically reasonable, so I'm removing them. 
# Begin by using filtered rmax dataset above.
my_152_site_list <- dat_out_rmed_Rhat$site_name

dat_out_cmed_rmax <- dat_out_cmed_2yrQ %>%
  filter(site_name %in% my_152_site_list) 

# And remove negative values.
dat_out_cmed_pos <- dat_out_cmed_rmax %>%
  filter(c_med > 0) # Removes 0 sites. Yay!

#### Rhat filter for c ####

# Before proceeding with the analyses, I will be filtering out sites at which
# the model did not converge well for the c parameter.
# Sites with Rhat > 1.05 will not pass muster.

dat_out_cmed_Rhat <- dat_out_cmed_pos %>%
  filter(Rhat < 1.05) # An additional 14 sites drop off.

# Export for future use.
saveRDS(dat_out_cmed_Rhat, "data_working/QcQ2_filtered_138sites_051023.rds")

# End of script.
