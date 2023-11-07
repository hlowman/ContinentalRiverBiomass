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

# The following script will calculate the maximum algal yields.
# Then, it will examine if they have any relationship with other covariates.

# Yield = (rmax*(0.5-(0.07*rmax))) / -lambda [based on Scheuerell 2016 PeerJ]

#### Setup ####

# Load packages.
lapply(c("lubridate","tidyverse", "here", "viridis",
         "reshape2","ggExtra","patchwork", "gsl"), require, character.only=T)

# Load datasets.

# Load in list containing all iterations of site-level parameters.
dat_out <- readRDS("data_working/beartooth_181rivers_model_params_all_iterations_050823.rds")

# And filtered rmax dataset for joining.
dat_rmax <- readRDS("data_working/rmax_filtered_143sites_110723.rds")

#### Format Data ####

# Take list containing all iterations of parameters and make into a df.
dat_out_df <- map_df(dat_out, ~as.data.frame(.x), .id="site_name")

#### Yield Formulas ####

# Equation to calculate max. algal yield is:
# yield = (r*(0.5-(0.07*r)))/-lambda

# Create new columns to calculate yield for EACH ITERATION at each site.
dat_out_df <- dat_out_df %>%
  mutate(yield = (r*(0.5-(0.07*r)))/-lambda)

# Note, I've added an extra negative value to the denominator in the second
# and third formulas since our lambda values are negative, but the manuscript
# assumes b values are positive (see Table 1, Scheuerell 2016 PeerJ).

# Calculate median yield values for each site and yield metric.
dat_out_yield_med <- dat_out_df %>%
  group_by(site_name) %>%
  summarize(yield_med = median(yield),
            yield_2.5 = quantile(yield, probs = 0.025),
            yield_97.5 = quantile(yield, probs = 0.975),
            lambda_med = median(lambda)) %>%
  ungroup() %>%
  # Some lower quantiles of yield were negative, so need to
  # filter those to become 0 instead so they plot.
  mutate(yield_2.5_ed = case_when(yield_2.5 < 0 ~ 0,
                                   TRUE ~ yield_2.5))

# Quick plot.
hist(dat_out_yield_med$yield_med) # still some negative values hmmm...
# Oh, this is because we don't yet have the poorly performing sites
# filtered out.

# Combine with rmax dataset.
# Note, this will trim down from 181 to 152 sites due to model diagnostics.
dat_out_yield_rmax <- left_join(dat_rmax, dat_out_yield_med)

# Export for future use.
saveRDS(dat_out_yield_rmax, "data_working/maxalgalyield_143sites_110723.rds")

# End of script.
