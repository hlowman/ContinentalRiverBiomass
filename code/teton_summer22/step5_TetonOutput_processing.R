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
# downloaded and added to a folder of the appropriate name prior to running the 
# code.

############################
## Setup
############################

## Load packages
lapply(c("tidyverse", "lubridate",
         "rstan","bayesplot","shinystan", "here"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################

# Import all five outputs from Teton. First four, because I'd split them
# to run in parallel and the fifth because on site misbehaved on it's 
# first run.
data_out1 <- readRDS("data_teton/teton_190rivers_output_pt1_2022_07_29.rds")
data_out2 <- readRDS("data_teton/teton_190rivers_output_pt2_2022_07_29.rds")
data_out3 <- readRDS("data_teton/teton_190rivers_output_pt3_2022_07_29.rds")
data_out4 <- readRDS("data_teton/teton_190rivers_output_pt4_2022_07_29.rds")
data_out5 <- readRDS("data_teton/teton_190rivers_output_pt1_1site_2022_08_10.rds")

# So, on the first attempt (8/16/2022), all of the above datasets pulled in empty.
# This was because there were some NA values in the light data, so I've gone back,
# and re-compiled/re-sent the data to Teton to fix this (8/17/2022).

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
data_out5_divs <- map(data_out5, extract_divergences)

# Make all of the above dataframes.
data_out1_divs_df <- map_df(data_out1_divs, 
                                   ~as.data.frame(.x), .id="site_name") %>%
  rename(divergences = `.x`)
data_out2_divs_df <- map_df(data_out2_divs, 
                                   ~as.data.frame(.x), .id="site_name") %>%
  rename(divergences = `.x`)
data_out3_divs_df <- map_df(data_out3_divs, 
                                   ~as.data.frame(.x), .id="site_name") %>%
  rename(divergences = `.x`)
data_out4_divs_df <- map_df(data_out4_divs, 
                                   ~as.data.frame(.x), .id="site_name") %>%
  rename(divergences = `.x`)
data_out5_divs_df <- map_df(data_out5_divs, 
                                   ~as.data.frame(.x), .id="site_name") %>%
  rename(divergences = `.x`)

# And join them together.
data_out_divs_all <- rbind(data_out1_divs_df, data_out2_divs_df)
data_out_divs_all <- rbind(data_out_divs_all, data_out3_divs_df)
data_out_divs_all <- rbind(data_out_divs_all, data_out4_divs_df)
data_out_divs_all <- rbind(data_out_divs_all, data_out5_divs_df)

# Export dataset
# saveRDS(data_out_divs_all, 
#         file = "data_working/teton_190rivers_model_divs_081622.rds")

# End of script.
