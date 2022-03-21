## Fitting models to data
## Step FOUR in Metabolism Modeling Workflow
## March 21, 2022
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to fit the Ricker model to a few sites of data to get the random effect working.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data - with new index number by site
df <- readRDS("data_working/df_207sites_indexed.rds")

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data
stan_data_compile <- function(x){
  data <- list(Ndays = x$Ndays, 
               light = x$light_rel, 
               GPP = x$GPP,
               GPP_sd = x$GPP_sd, 
               tQ = x$tQ,
               new_e = x$new_e) # new column for reinitialization of B[j] values
  
  return(data)
}

stan_data_l <- lapply(df, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model

# sets initial values of c and s to help chains converge
# init_Ricker <- function(...) {
#   list(c = 0.5, s = 0.5) # new values as of jan 2022
# }

## export results

PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("/project/modelscape/users/hlowman/jobscripts/teton_moresites/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_noP.stan",
                                                data = x, chains = 3,iter = 5000,
                                                #init = init_Ricker,
                                                control = list(max_treedepth = 12)))

saveRDS(PM_outputlist_Ricker, "/project/modelscape/users/hlowman/jobresults/teton_moresites/teton_131rivers_output_Ricker_2022_02_07.rds")

# End of script.