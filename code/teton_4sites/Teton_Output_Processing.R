## 4 river Teton output script
## January 21, 2022
## Heili Lowman

# The following script will briefly examine the output of the
# small run I sent to Teton with the new edits made to the model
# structure, namely (1) reinitialization following gaps in time
# series > 14 days and (2) choosing between models with and without
# P terms in an effort to improve fit.

# Load packages
lapply(c("calecopal", "cowplot", "parallel",
         "devtools", "viridis", "here",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", 
         "ggrepel", "patchwork", "MCMCglmm",
         "data.table", "rlist", "pipeR"), require, character.only=T)

#### Data Import & Processing ####

teton4 <- readRDS("data_teton/teton_4rivers_output_Ricker_2022_01_22.rds")

#### Extraction of model parameters ####

# Extract the parameters resulting from fitting the data to the model.

# Going to create a function of the above to pull out data of interest from
# all sites.
extract_params <- function(df){
  extract(df, c("r","lambda","s","c","sig_p","sig_o"))
}

# And now map this to the entire output list.
data_out_params <- map(teton4, extract_params)
# P, s, and c should not exist for 2 sites.
# Yay! They don't exist!! a.k.a. I get an error message.

# So, I actually will need to modify the extract function based on what
# parameters do exist. This will need to be based on a better filter in the
# future, but for now, I'm going to hard-wire it.

# with P
data_out_params_wP <- teton4[c("nwis_07191222", "nwis_04121944")] %>>%
  map(extract_params)
  
# create new function for extracting parameters without P present
extract_params_woP <- function(df){
  extract(df, c("r","lambda","sig_p","sig_o"))
}

# without P
data_out_params_woP <- teton4[c("nwis_02266300", "nwis_03067510")] %>>%
  map(extract_params_woP)

# And create a dataframe
params_df1 <- map_df(data_out_params_wP, ~as.data.frame(.x), .id="site_name") %>%
  # and add "K" to it, calculating for each individual iteration
  mutate(k = (-1*r)/lambda)

params_df2 <- map_df(data_out_params_woP, ~as.data.frame(.x), .id="site_name") %>%
  # and add "K" to it, calculating for each individual iteration
  mutate(k = (-1*r)/lambda) %>%
  mutate(c = NA,
         s = NA) %>%
  select(site_name, r, lambda, s, c, sig_p, sig_o, k)

params_join <- bind_rows(params_df1, params_df2)

# And now to calculate means by site.
sim_sums <- params_join %>%
  group_by(site_name) %>%
  summarize(r_mean = mean(r),
            k_mean = mean(k),
            lambda_mean = mean(lambda),
            s_mean = mean(s),
            c_mean = mean(c),
            sigp_mean = mean(sig_p),
            sigo_mean = mean(sig_o)) %>%
  ungroup()

launch_shinystan(teton4$nwis_02266300) # No divergences!!
launch_shinystan(teton4$nwis_03067510) # Very few divergences!! (4)
launch_shinystan(teton4$nwis_04121944) # Same! (18)
launch_shinystan(teton4$nwis_07191222) # Same! (7)

# End of script.
