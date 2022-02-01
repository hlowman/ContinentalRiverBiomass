## Fitting models to data - JR Blaszczak
## Modified: June 23, 2021
## Heili Lowman

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here", "viridis"), require, character.only=T)

## Source data
df <- readRDS("NWIS_1site_subset_SL.rds")

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
               tQ = x$tQ)
  return(data)
}

stan_data_l <- stan_data_compile(NWIS_1site_subset_SL)

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

## PM 2 - Latent Biomass (Ricker)
# With Persistence Term (P)
init_Ricker <- function(...) {
  list(c = 0.5, s = 0.5)
}

## export results
PM_outputlist_Ricker <- stan("Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                            data = stan_data_l,
                            chains = 3,iter = 5000,
                            init = init_Ricker,
                            control = list(max_treedepth = 12))

saveRDS(PM_outputlist_Ricker, "stan_1riv_output_Ricker_2021_06_23.rds")

# End of script.
