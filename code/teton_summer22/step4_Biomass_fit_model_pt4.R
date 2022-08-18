## Resilience of Stream Productivity to Disturbance
## July 29, 2022
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

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data
df <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_summer22/list_190sites_10yrQnorm_allSL_pt4.rds")

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile necessary data
stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), # number of days/records
               light = x$light_rel, # relativized measured stream light
               GPP = x$GPP,         # GPP estimates
               GPP_sd = x$GPP_sd,   # standard deviation of GPP estimates
               tQ = x$Q_rel,        # relativized discharge (to 10yr flood)
               new_e = x$new_e)     # instances where re-initialization needed
               
  return(data)
}

stan_data_l <- lapply(df, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model

# sets initial values to help chains converge
init_Ricker <- function(...) {
  list(r = 0.2, lambda = -0.03, c = 0.5, s = 1.5) # values to match priors
}

## export results
PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("/project/modelscape/users/hlowman/jobscripts/teton_summer22/Stan_ProductivityModel.stan",
                                                data = x, chains = 3,iter = 5000,
                                                init = init_Ricker,
                                                control = list(max_treedepth = 12)))

saveRDS(PM_outputlist_Ricker, "/project/modelscape/users/hlowman/jobresults/teton_summer22/teton_190rivers_output_pt4_2022_08_18.rds")

# End of script.
