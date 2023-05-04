## Resilience of Stream Productivity to Disturbance
## Originally created: October 12, 2022
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

# Please also note, all file directories below are for the Beartooth HPC hosted
# by the University of Wyoming ARCC. All models are run there with outputs then
# transferred back to this working directory.

#### Setup ####

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data
df <- readRDS("/project/modelscape/users/hlowman/jobscripts/beartooth_spring23/list_181sites_Qmaxnorm_SavoySL_pt2.rds")

#### Stan data prep ####

rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile necessary data
stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), # number of days/records
               light = x$light_rel, # relativized light (to maximum light)
               GPP = x$GPP,         # GPP estimates
               GPP_sd = x$GPP_sd,   # standard deviation of GPP estimates
               tQ = x$Q_rel,        # relativized discharge (to maximum Q)
               new_e = x$new_e)     # instances where re-initialization needed
               
  return(data)
}

stan_data_l <- lapply(df, function(x) stan_data_compile(x))

#### Run Stan to get parameter estimates ####

# Latent Biomass (Ricker population) Model

# sets initial values to help chains converge
init_Ricker <- function(...) {
  list(r = 0.2, lambda = -0.03, c = 0.5, s = 1.5) # values to match priors
}

## export results
PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("/project/modelscape/users/hlowman/jobscripts/beartooth_spring23/Stan_ProductivityModel.stan",
                                                data = x, 
                                                chains = 3,
                                                iter = 5000,
                                                init = init_Ricker,
                                                control = list(max_treedepth = 12)))

saveRDS(PM_outputlist_Ricker, "/project/modelscape/users/hlowman/jobresults/beartooth23/teton_181rivers_output_pt2_2032_05_04.rds")

# End of script.
