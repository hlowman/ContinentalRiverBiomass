## Fitting models to data
## Created: October 13, 2021
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to fit the Ricker model to 3 sites of data.

# I'm going to run this same code twice - first to see the times series
# as they currently are and second to see how they do with gaps
# (a.k.a. data I've purposefully removed).

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here", "viridis"), 
       require, character.only=T)

#### Run 1 ####

## Source data
df <- readRDS("data_working/df_3sites.rds")

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

stan_data_l <- lapply(df, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model
# With Persistence Term (P)

# sets initial values of c and s to help chains converge
init_Ricker <- function(...) {
  list(c = 0.5, s = 100)
}

## export results
PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("code/gap_experimenting/Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                                                data = x,chains = 3,iter = 5000,
                                                init = init_Ricker,
                                                control = list(max_treedepth = 12)))

saveRDS(PM_outputlist_Ricker, "data_working/stan_3rivers_output_Ricker_2021_10_13.rds")

#### Run 2 ####
# some of the functions from above will be re-used
## Source data
df_g <- readRDS("data_working/df_3sites_gappy.rds")

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

stan_data_l_g <- lapply(df_g, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model
# With Persistence Term (P)

## export results
PM_outputlist_Ricker_g <- lapply(stan_data_l_g,
                               function(x) stan("code/gap_experimenting/Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                                                data = x,chains = 3,iter = 5000,
                                                init = init_Ricker,
                                                control = list(max_treedepth = 12)))

saveRDS(PM_outputlist_Ricker_g, "data_working/stan_3rivers_gaps_output_Ricker_2021_10_13.rds")

#### Run 3 ####
# again repurposing some of the functions from above
## Source data
df_y <- readRDS("data_working/df_1site_allyrs.rds")

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

stan_data_l_y <- lapply(df_y, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model
# With Persistence Term (P)

## export results
PM_outputlist_Ricker_y <- lapply(stan_data_l_y,
                                 function(x) stan("code/gap_experimenting/Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                                                  data = x,chains = 3,iter = 5000,
                                                  init = init_Ricker,
                                                  control = list(max_treedepth = 12)))

saveRDS(PM_outputlist_Ricker, "data_working/stan_1river_years_output_Ricker_2021_10_13.rds")

# End of script.
