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

# This data was part of a re-run where Q was normalized to Qmax instead of Q10yr flood.
data_out5 <- readRDS("data_teton/teton_4rivers_output_2022_09_08.rds")
data_out6 <- readRDS("data_teton/teton_2rivers_output_2022_09_12.rds")

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

# Note - I crashed the desktop trying to run this, so only do the below
# using the Pinyon server !!!

# Going to create a function to pull out all iterations.

extract_all <- function(df){
  extract(df, c("r", "lambda", "s", "c", 
                "B", "P", "pred_GPP", "sig_p", "sig_o"))
}

# And now map this to all output lists. 
# (Be patient - this step takes a while.)
data_out1_its <- map(data_out1, extract_all)
data_out2_its <- map(data_out2, extract_all)
data_out3_its <- map(data_out3, extract_all)
data_out4_its <- map(data_out4, extract_all)

# And for re-runs
data_out5_its <- map(data_out5, extract_all)
data_out6_its <- map(data_out6, extract_all)

# Saving these four out just in case because the load-in took so long.
# saveRDS(data_out1_its,
#        file = "data_working/teton_190rivers_model_all_params_all_iterations_082322pt1.rds")
# saveRDS(data_out2_its,
#        file = "data_working/teton_190rivers_model_all_params_all_iterations_082322pt2.rds")
# saveRDS(data_out3_its,
#        file = "data_working/teton_190rivers_model_all_params_all_iterations_082322pt3.rds")
# saveRDS(data_out4_its,
#        file = "data_working/teton_190rivers_model_all_params_all_iterations_082322pt4.rds")

# So, these can't be bound together because I've extracted "B" and "pred_GPP"
# which have different numbers of days and prevent the rbind() function below
# from working.

# Instead, I'm going to proceed with extracting only the "r" values for now.

#### Load-in for Parameter Extraction ####

# Load in the datasets saved above.
data_out1_its <- readRDS("data_working/teton_190rivers_model_all_params_all_iterations_082322pt1.rds")
data_out2_its <- readRDS("data_working/teton_190rivers_model_all_params_all_iterations_082322pt2.rds")
data_out3_its <- readRDS("data_working/teton_190rivers_model_all_params_all_iterations_082322pt3.rds")
data_out4_its <- readRDS("data_working/teton_190rivers_model_all_params_all_iterations_082322pt4.rds")

#### r, lambda, c, s Parameters ####

# Select the site-level parameter posterior predictions only.
# create a function to do so.
param_compile <- function(x){
  data <- list(r = x$r,         # maximum growth rate
               lambda = x$lambda,   # carrying capacity
               c = x$c,        # critical discharge to scour biomass
               s = x$s)     # sensitivity of persistence curve
  
  return(data)
}

# And apply to each of the four lists.
data_out1_its_params <- lapply(data_out1_its, function(x) param_compile(x))
data_out2_its_params <- lapply(data_out2_its, function(x) param_compile(x))
data_out3_its_params <- lapply(data_out3_its, function(x) param_compile(x))
data_out4_its_params <- lapply(data_out4_its, function(x) param_compile(x))

# And for re-run.
data_out5_its_params <- lapply(data_out5_its, function(x) param_compile(x))
data_out6_its_params <- lapply(data_out6_its, function(x) param_compile(x))

# And create dataframes.
data_out1_its_pdf <- map_df(data_out1_its_params, ~as.data.frame(.x), .id="site_name")
data_out2_its_pdf <- map_df(data_out2_its_params, ~as.data.frame(.x), .id="site_name")
data_out3_its_pdf <- map_df(data_out3_its_params, ~as.data.frame(.x), .id="site_name")
data_out4_its_pdf <- map_df(data_out4_its_params, ~as.data.frame(.x), .id="site_name")

# And join them together.
data_out_its_pall <- rbind(data_out1_its_pdf, data_out2_its_pdf)
data_out_its_pall <- rbind(data_out_its_pall, data_out3_its_pdf)
data_out_its_pall <- rbind(data_out_its_pall, data_out4_its_pdf)

# And for re-run.
data_out_its_4sites <- map_df(data_out5_its_params, ~as.data.frame(.x), .id="site_name")
data_out_its_2sites <- map_df(data_out6_its_params, ~as.data.frame(.x), .id="site_name")

# the above line of code sometimes doesn't play nicely if R has been running
# for awhile, so the fix is to exit RStudio and reopen the project/file
# OR, as is the case with this, where data_out has just taken an hour to load,
# you should instead uncheck and re-check 'rstan' so that extract function
# takes precedence over the same function in the 'tidyr' package.

# Export data.
# saveRDS(data_out_its_pall,
#        file = "data_working/teton_190rivers_model_site_params_all_iterations_082422.rds")

# saveRDS(data_out_its_4sites,
#        file = "data_working/teton_4rivers_model_site_params_all_iterations_090822.rds")

# saveRDS(data_out_its_2sites,
#        file = "data_working/teton_2rivers_model_site_params_all_iterations_091222.rds")

# End of script.
