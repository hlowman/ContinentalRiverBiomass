## Fitting models to data
## January 27, 2022
## Heili Lowman

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data
df <- readRDS("NWIS_1site_subset_gaps.rds")

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data
stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), 
               light = x$light_rel, 
               GPP = x$GPP,
               GPP_sd = x$GPP_sd, 
               tQ = x$tQ,
               new_e = x$new_e) # new column for reinitialization of B[j] values
               
  return(data)
}

stan_data_l <- stan_data_compile(df)

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model

# sets initial values of c and s to help chains converge
init_Ricker <- function(...) {
  list(c = 0.5, s = 0.5) # new values as of jan 2022
}

## export results

PM_outputlist_Ricker <- stan("Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP.stan",
                             data = stan_data_l,
                             chains = 3,iter = 5000,
                             init = init_Ricker,
                             control = list(max_treedepth = 12))

saveRDS(PM_outputlist_Ricker, "output_Ricker_2022_02_23.rds")

#### BONUS CODE ####

########################################
# for re-initialization column "new_e" #
########################################

# Addition of time-series delineation column using seqFUN function created below.
# NOTE - this was done prior to the creation of the dataframe "df" used above.

# create placeholder columns for day to day differences and sequences
data <- data %>%
  mutate(diff_time = 0, seq = 1, new_e = 0)

## split list by ID
l <- split(data, data$site_name)

# loop over the separate time sequences for a given site
seqFUN <- function(d){
  
  # calculate the difference from one day to the next
  for(i in 2:nrow(d)){
    d$diff_time[i] = difftime(time1 = d$date[i], time2 = d$date[(i-1)],
                              units = "days")
  }
  
  # delineate sequenced time frames based on day to day differences
  # anything less than 14 day gaps is permissible
  for(i in 2:nrow(d)){
    if(d$diff_time[i] < 14){
      d$seq[i] = d$seq[(i-1)]
    } else {
      d$seq[i] = d$seq[(i-1)] + 1
    }
  }
  
  # add column to delineate changes in events
  for(i in 2:nrow(d)){
    d$new_e[i] = d$seq[i] - d$seq[i-1]
  }
  
  return(d)
  
}

# apply event delineation function
l <- lapply(l, function(x) seqFUN(x))

####################################################
# extracting parameters/divergences/diagnostics... #
####################################################

# After fitting the model above, here is code to extract information of interest.
data_out <- readRDS("output_Ricker_2022_02_23.rds")

#### Parameters ####

# Going to create a function to pull out data of interest.
# Sticking to just the four parameters of interest for initial data viz.
extract_params <- function(df){
  extract(df, c("r","lambda","s","c"))
}

data_out_params <- extract_params(data_out)
# the above line of code sometimes doesn't play nicely if R has been up and running
# for awhile, so the fix is to exit RStudio and reopen the project/file
# OR, you should instead uncheck and re-check 'rstan' so that it's extract function
# takes precedence over the same function in the 'tidyr' package.

# And add "K" to it, calculating for each individual iteration
data_out_params_df <- data_out_params %>%
  mutate(k = (-1*r)/lambda)

#### Divergences ####

# Also, need to extract divergence information from the sites.
# A method for extrating divergences per the stan handbook:
# mc-stan.org/rstan/reference/check_hmc_diagnostics.html
extract_divergences <- function(df){
  divergences <- as.numeric(get_num_divergent(df))
}

# this function can also be finicky, and require you to quit R and re-load everything back in.
data_out_divergences <- extract_divergences(data_out)
# WOOHOO!!

#### Check model diagnostics ####

# Using rstan documentation found at:
# https://mc-stan.org/rstan/articles/stanfit_objects.html

# In $summary, results for all chains are merged
# se_mean = Monte Carlo standard error
# n_eff = effective sample size
# Rhat = R-hat statistic

# Going to create a function to pull out summaries from all sites
extract_summary <- function(x){
  df <- x
  df1 <- summary(df,
                 pars = c("r", "lambda", "s", "c"))$summary
  as.data.frame(df1) %>% rownames_to_column("parameter")
}

# (Be patient - This step takes a while.)
data_out_summary <- extract_summary(data_out)

# End of script.
