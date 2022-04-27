## 206 rivers data output from Teton
## Created: April 22, 2022
## Heili Lowman

# The following code will do some preliminary processing of the output
# of the 206 site dataset sent to Teton on April 22, 2022.

# The visualization code has been moved to TetonOutput_exploration.R.

# Load packages
lapply(c("calecopal", "cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", "here",
         "ggrepel", "patchwork", "grid","gridExtra"), require, character.only=T)

#### Data Import & Processing ####

data_out <- readRDS("data_teton/teton_206rivers_output_Ricker_2022_04_22_2.rds")

# Make sure we can extract parameters by testing at one site.
test_params <- extract(data_out$nwis_01608500, c("r","lambda","s","c"))
# Yay - this works!

# Going to create a function of the above to pull out data of interest from
# all sites. Sticking to just the four parameters of interest for initial data viz.
extract_params <- function(df){
  extract(df, c("r","lambda","s","c"))
}

# And now map this to the entire output list.
data_out_params <- map(data_out, extract_params)
# the above line of code sometimes doesn't play nicely if R has been up and running
# for awhile, so the fix is to exit RStudio and reopen the project/file
# OR, as is the case with this code, where data_out has just taken an hour to load,
# you should instead uncheck and re-check 'rstan' so that it's extract function
# takes precedence over the same function in the 'tidyr' package.

# And create a dataframe
data_out_params_df <- map_df(data_out_params, ~as.data.frame(.x), .id="site_name") %>%
  # and add "K" to it, calculating for each individual iteration
  mutate(k = (-1*r)/lambda)

# Export dataset
# saveRDS(data_out_params_df, 
#        file = "data_working/teton_207rivers_model_parameters_all_iterations_020622.rds")