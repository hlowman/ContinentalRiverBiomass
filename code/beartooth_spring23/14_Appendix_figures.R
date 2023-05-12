## Recovery of Stream Productivity following Disturbance
## Originally created: May 12, 2023
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
lapply(c("tidybayes", "brms", "tidyverse", "lubridate", 
         "data.table", "GGally", "patchwork", "here"), 
       require, character.only=T)

#### Data ####

# Import necessary datasets.

# First, the data for the maximum accrual (amax) models.
dat_amax <- readRDS("data_working/amax_covariates_152sites_051123.rds")

# Next, the data for the Qc:Q2yr models.
dat_Qc <- readRDS("data_working/Qc_covariates_138sites_051123.rds")

# And, the original data used in modeling.
dat_in <- readRDS("data_working/df_181sites_Qmaxnorm_SavoySL.rds")

#### Timeseries Length Appendix Figure ####

# Making a histogram of the lengths of the timeseries for an appendix figure.
site_lengths <- dat_in %>%
  group_by(site_name) %>%
  summarize(days = n()) %>% # count the number of rows/days
  ungroup()

# Adding bins for more intuitive plotting.
site_lengths <- site_lengths %>%
  mutate(bin = factor(case_when(days < 180 ~ "3-6 months",
                                days >= 180 & days < 365 ~ "6 mos. - 1 year",
                                days >= 365 & days < 730 ~ "1-2 years",
                                days >= 730 & days < 1460 ~ "2-4 years",
                                days >= 1460 & days < 2190 ~ "4-6 years",
                                days >= 2190 & days < 2920 ~ "6-8 years",
                                days >= 2920 ~ "8-10 years"),
                      levels = c("3-6 months", "6 mos. - 1 year",
                                 "1-2 years", "2-4 years",
                                 "4-6 years", "6-8 years",
                                 "8-10 years")))

(ts_hist <- ggplot(site_lengths, aes(x = bin)) +
    geom_histogram(stat = "count", fill = "#7AC9B7") +
    labs(x = "Length of Timeseries",
         y = "Number of Sites") +
    theme_bw() +
    theme(text = element_text(size = 10)))

# And export.
# ggsave(ts_hist,
#        filename = "figures/beartooth_spring23/TS_length_hist_051223.jpg",
#        width = 14,
#        height = 7,
#        units = "cm")

# And calculate median timeseries length for inclusion in methods.
median(site_lengths$days) # 870

# End of script.

