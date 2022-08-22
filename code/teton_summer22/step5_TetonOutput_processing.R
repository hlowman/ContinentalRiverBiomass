## Resilience of Stream Productivity to Disturbance
## August 16, 2022
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

#### Setup ####

## Load packages
lapply(c("tidyverse", "lubridate",
         "rstan","bayesplot","shinystan", "here"), require, character.only=T)

#### Data Import ####

# Import all outputs from Teton.
data_out1 <- readRDS("data_teton/teton_190rivers_output_pt1_2022_08_17.rds")
data_out2 <- readRDS("data_teton/teton_190rivers_output_pt2_2022_08_18.rds")
data_out3 <- readRDS("data_teton/teton_190rivers_output_pt3_2022_08_18.rds")
data_out4 <- readRDS("data_teton/teton_190rivers_output_pt4_2022_08_18.rds")

# So, on the first attempt (8/16/2022), the above datasets pulled in empty.
# This was because there were some NA values in the light data, 
# so I've gone back and re-compiled/re-sent the data to Teton (8/17/2022).

#### Divergences ####

# First, need to pull out divergences from all datasets.
# so, creating a function to do so...
extract_divergences <- function(df){
  divergences <- as.numeric(get_num_divergent(df))
}

# ...and now map this to the output datasets.
# This function may be finicky and require you to quit R and re-load.
data_out1_divs <- map(data_out1, extract_divergences)
data_out2_divs <- map(data_out2, extract_divergences)
data_out3_divs <- map(data_out3, extract_divergences)
data_out4_divs <- map(data_out4, extract_divergences)

# Make all of the above dataframes.
data_out1_divs_df <- map_df(data_out1_divs, 
                                   ~as.data.frame(.x), .id="site_name") %>%
  rename(divergences = `.x`) %>%
  mutate(dataset = 1)
data_out2_divs_df <- map_df(data_out2_divs, 
                                   ~as.data.frame(.x), .id="site_name") %>%
  rename(divergences = `.x`) %>%
  mutate(dataset = 2)
data_out3_divs_df <- map_df(data_out3_divs, 
                                   ~as.data.frame(.x), .id="site_name") %>%
  rename(divergences = `.x`) %>%
  mutate(dataset = 3)
data_out4_divs_df <- map_df(data_out4_divs, 
                                   ~as.data.frame(.x), .id="site_name") %>%
  rename(divergences = `.x`) %>%
  mutate(dataset = 4)

# And join them together.
data_out_divs_all <- rbind(data_out1_divs_df, data_out2_divs_df)
data_out_divs_all <- rbind(data_out_divs_all, data_out3_divs_df)
data_out_divs_all <- rbind(data_out_divs_all, data_out4_divs_df)

# Export dataset
# saveRDS(data_out_divs_all,
#         file = "data_working/teton_190rivers_model_divs_081922.rds")

# 4 above 1000 - removing these because upon closer shiny stan inspection, n_eff values also tended to be quite low
# 4 above 500 - these looked alright except the 0112400 site
# 8 above 300
# 17 above 100 - checked a few sporadically and Rhat and n_eff values were all over the place, so I need to examine model diagnostics more closely it seems

#### Diagnostics ####

# Using rstan documentation found at:
# https://mc-stan.org/rstan/articles/stanfit_objects.html

# Obtain summary statistics using new "extract_summary2" function.
extract_summary <- function(x){
  df <- x
  df1 <- summary(df,
                 pars = c("r", "lambda", "s", "c"),
                 probs = c(0.025, 0.975))$summary # 2.5% and 97.5% percentiles as bounds for the 95% confidence interval, as reported by Appling
  as.data.frame(df1) %>% rownames_to_column("parameter")
}

# And now map this to the entire output list. 
# (Be patient - this step takes a while.)
data_out1_diags <- map(data_out1, extract_summary)
data_out2_diags <- map(data_out2, extract_summary)
data_out3_diags <- map(data_out3, extract_summary)
data_out4_diags <- map(data_out4, extract_summary)

# And create a dataframe
data_out1_diags_df <- map_df(data_out1_diags, 
                               ~as.data.frame(.x), .id="site_name")
data_out2_diags_df <- map_df(data_out2_diags, 
                               ~as.data.frame(.x), .id="site_name")
data_out3_diags_df <- map_df(data_out3_diags, 
                               ~as.data.frame(.x), .id="site_name")
data_out4_diags_df <- map_df(data_out4_diags, 
                               ~as.data.frame(.x), .id="site_name")

# And join them together.
data_out_diags_all <- rbind(data_out1_diags_df, data_out2_diags_df)
data_out_diags_all <- rbind(data_out_diags_all, data_out3_diags_df)
data_out_diags_all <- rbind(data_out_diags_all, data_out4_diags_df)

# Export dataset
# saveRDS(data_out_diags_all,
#         file = "data_working/teton_190rivers_model_diags_081922.rds")

# And now to institute filtering to select sites for final analysis.
# For parameters "r", "lambda", "c", and "s",
# - Rhat < 1.05
# - n_eff > 750 (10% of total 7,500 non-warm-up iterations)
# For full model,
# - # of divergences < 100

data_out_filter <- data_out_diags_all %>%
  filter(Rhat < 1.05) %>% # 5 sites drop off
  filter(n_eff > 750) # 16 sites drop off

# Going to take a closer look to see just how many sites have these red
# flags for various parameters

data_out_out1 <- data_out_diags_all %>%
  filter(Rhat >= 1.05) # c appears to be causing the most issues

sites_out_Rhat <- unique(data_out_out1$site_name)

data_out_out2 <- data_out_diags_all %>%
  filter(n_eff <= 750) # c again appears to be causing the most issues

sites_out_neff <- unique(data_out_out2$site_name)

data_out_out3 <- data_out_divs_all %>%
  filter(divergences > 100)

sites_out_divs <- unique(data_out_out3$site_name)

# And create a dataset that compiles these summary statistics.
data_sites <- data_out_divs_all %>%
  # Create a column with full listing of 190 sites
  select(site_name) %>%
  # Add column with "Y" for sites below 100 divergences, "N" for sites above
  mutate(divs_less_100 = if_else(site_name %in% sites_out_divs, "N", "Y")) %>%
  # Add a column with "Y" for sites with all Rhat < 1.05, "N" for sites above
  mutate(Rhat_less_1.05 = if_else(site_name %in% sites_out_Rhat, "N", "Y")) %>%
  # Add column with "Y" for sites with all n_eff > 750, "N" for sites below
  mutate(neff_more_750 = if_else(site_name %in% sites_out_neff, "N", "Y"))

ally <- data_sites %>%
  filter(divs_less_100 == "Y") %>%
  filter(Rhat_less_1.05 == "Y") %>%
  filter(neff_more_750 == "Y") # 97 sites

alln <- data_sites %>%
  filter(divs_less_100 == "N") %>%
  filter(Rhat_less_1.05 == "N") %>%
  filter(neff_more_750 == "N") # 10 sites

# For the time being, I am going to remove the 10 sites that pass none
# of the three rules.

sites_to_remove <- unique(alln$site_name)

# [1] "nwis_0165389205" "nwis_02171645"   "nwis_02336300"   "nwis_02336728"  
# [5] "nwis_03067510"   "nwis_03073000"   "nwis_03289200"   "nwis_03302030"  
# [9] "nwis_04124000"   "nwis_06893350" 

#### Parameter Posteriors ####



# End of script.
