## Fitting models to data
## July 13, 2021
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to practice fitting the Ricker model to ~30 site-years of data.
# The compiled data will then be run by sending the job to Teton (UWyoming).

# I've commented out those steps that I feel, for the time being, I don't
# need to perform, and I've changed the appropriate filepaths to match my
# repository structure.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here", "viridis"), 
       require, character.only=T)

## Source data - sources the file itself
source("code/teton_34sites/DataSource_34rivers_StreamLight.R")

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)#parallel::detectCores())

## compile data
stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), 
               light = x$light_rel, 
               GPP = x$GPP,
               GPP_sd = x$GPP_sd, 
               tQ = x$tQ)
  return(data)
}

stan_data_l <- lapply(df, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

## PM 2 - Latent Biomass (Ricker)
# With Persistence Term (P)
# sets initial values of c and s to help chain converge
init_Ricker <- function(...) {
  list(c = 0.5, s = 100)
}

## export results
PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                                                data=x,chains=3,iter=5000,
                                                init = init_Ricker,
                                                control=list(max_treedepth=12)))


# not going to save the elapsed time for Teton runs
#PM_Ricker_elapsedtime <- lapply(PM_outputlist_Ricker, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_Ricker, "/project/modelscape/users/hlowman/jobresults/stan_34rivers_output_Ricker_2021_07_14.rds")
#saveRDS(PM_Ricker_elapsedtime, "data_working/stan_34rivers_Ricker_time_2021_07_14.rds")

# End of script.
