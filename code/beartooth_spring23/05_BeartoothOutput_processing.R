## Recovery of Stream Productivity following Disturbance
## Originally created: October 13, 2022
## Heili Lowman

#### READ ME ####

# The following set of scripts will walk through the steps necessary to
# prep and send data to Beartooth as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the data loaded in this script is too large to be posted to GitHub.
# Instead, the outputs of this script are those that will be published with the
# manuscript at this link: INSERT ZENODO LINK HERE.

# Note - I CRASHED the desktop trying to run this, so only do the below
# using the Pinyon SERVER!!!

#### Setup ####

## Load packages
lapply(c("tidyverse", "lubridate",
         "rstan", "bayesplot",
         "shinystan", "here"), require, character.only=T)

#### Data Import ####

# Import all outputs from Beartooth.
data_out1 <- readRDS("data_beartooth/beartooth_181rivers_output_pt1_2023_05_05.rds")
data_out2 <- readRDS("data_beartooth/beartooth_181rivers_output_pt2_2023_05_05.rds")
data_out3 <- readRDS("data_beartooth/beartooth_181rivers_output_pt3_2023_05_05.rds")
data_out4 <- readRDS("data_beartooth/beartooth_181rivers_output_pt4_2023_05_05.rds")

#### Parameter Posteriors ####

##### Daily Estimates #####

# First, pull out all iterations of all parameters.

# Going to create a function to pull out all iterations.
extract_allparams <- function(df){
  rstan::extract(df, c("r", "lambda", "s", "c", 
                "B", "P", "pred_GPP", "sig_p", "sig_o"))
}

# And now map this to output lists. 
# (Be patient - this step takes a while.)
data_out1_its <- map(data_out1, extract_allparams)
data_out2_its <- map(data_out2, extract_allparams)
data_out3_its <- map(data_out3, extract_allparams)
data_out4_its <- map(data_out4, extract_allparams)

# Saving these out because the load-in took so long.
# saveRDS(data_out1_its,
#        file = "data_working/beartooth_181rivers_model_all_params_all_iterations_050823pt1.rds")
# saveRDS(data_out2_its,
#        file = "data_working/beartooth_181rivers_model_all_params_all_iterations_050823pt2.rds")
# saveRDS(data_out3_its,
#        file = "data_working/beartooth_181rivers_model_all_params_all_iterations_050823pt3.rds")
# saveRDS(data_out4_its,
#        file = "data_working/beartooth_181rivers_model_all_params_all_iterations_050823pt4.rds")

# These can't be bound together because I've extracted "B" and "pred_GPP"
# which have different numbers of days and prevent the rbind() function below
# from working.

##### Site-Level Estimates #####

# So, I'm going to proceed with extracting only the r, lambda, s, c, sig_p,
# sig_o values - the site-level parameters, since these are needed for the re-
# prediction of GPP to calculate RMSE values (and the datasets are smaller).

# Going to create a function to pull out all iterations.
extract_siteparams <- function(df){
  rstan::extract(df, c("r", "lambda", "s", "c", "sig_p", "sig_o"))
}

# And now map this to output lists.
data_out1_its_params <- map(data_out1, extract_siteparams)
data_out2_its_params <- map(data_out2, extract_siteparams)
data_out3_its_params <- map(data_out3, extract_siteparams)
data_out4_its_params <- map(data_out4, extract_siteparams)

# Concatenate lists.
data_out_its_params_all <- c(data_out1_its_params, data_out2_its_params, 
                             data_out3_its_params, data_out4_its_params)

# Saving this out too, but keeping as a list for iteration purposes.
# saveRDS(data_out_its_params_all,
#        file = "data_working/beartooth_181rivers_model_params_all_iterations_050823.rds")

#### Diagnostics ####

# Using rstan documentation found at:
# https://mc-stan.org/rstan/articles/stanfit_objects.html

# Obtain summary statistics using new "extract_summary" function.
# Note, this pulls diagnostics for all four site-level parameters.
extract_summary <- function(x){
  df <- x
  df1 <- summary(df,
                 pars = c("r", "lambda", "s", "c"),
                 probs = c(0.025, 0.975))$summary # 2.5% and 97.5% percentiles
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
data_out_diags_12 <- rbind(data_out1_diags_df, data_out2_diags_df)
data_out_diags_123 <- rbind(data_out_diags_12, data_out3_diags_df)
data_out_diags_all <- rbind(data_out_diags_123, data_out4_diags_df)

# Export dataset
# saveRDS(data_out_diags_all,
#         file = "data_working/beartooth_181rivers_model_diags_050923.rds")

#### Divergences ####

# Create a function to extract divergences.
extract_divergences <- function(df){
  divergences <- as.numeric(get_num_divergent(df))
}

# ...and now map this to the output datasets.
# This function may be finicky and require you to quit R and re-load.
# So, be sure you run this last.
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
  mutate(dataset = 1)
data_out3_divs_df <- map_df(data_out3_divs, 
                            ~as.data.frame(.x), .id="site_name") %>%
  rename(divergences = `.x`) %>%
  mutate(dataset = 1)
data_out4_divs_df <- map_df(data_out4_divs, 
                            ~as.data.frame(.x), .id="site_name") %>%
  rename(divergences = `.x`) %>%
  mutate(dataset = 1)

# And join them together.
data_out_divs_12 <- rbind(data_out1_divs_df, data_out2_divs_df)
data_out_divs_123 <- rbind(data_out_divs_12, data_out3_divs_df)
data_out_divs_all <- rbind(data_out_divs_123, data_out4_divs_df)

# Export dataset
# saveRDS(data_out_divs_all,
#         file = "data_working/beartooth_181rivers_model_divs_050923.rds")

# End of script.
