## Fitting models to data
## Step FOUR in Metabolism Modeling Workflow
## February 13, 2022
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to fit the Ricker model to 1 site of data without a P term. This
# site was for some reason breaking the Teton run, so for time's sake,
# I'm running it separately here.

# I double-checked all the code steps prior to this one on 01/27/22.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data
df <- readRDS("data_working/df_207sites.rds")

# Pull out only site of interest
df1 <- df$nwis_01201487

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

stan_data_l <- stan_data_compile(df1)

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model

# sets initial values of c and s to help chains converge
# init_Ricker <- function(...) {
#   list(c = 0.5, s = 0.5) # new values as of jan 2022
# }

## export results

PM_outputlist_Ricker <- stan(here("code", "teton_moresites", "Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_noP.stan"),
                                                data = stan_data_l, chains = 3,iter = 5000,
                                                control = list(max_treedepth = 12))

# STILL GETTING AN ERROR MESSAGE, SO, THE SITE WILL BE PRESENTED ONLY FROM THE RUN INCLUDING THE P TERM.
# nwis_01201487 = Still River, Brookfield Center, CT

# End of script.
