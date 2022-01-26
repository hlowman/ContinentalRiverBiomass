## Fitting models to data
## Step THREE in Metabolism Modeling Workflow
## January 26, 2022
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to fit the Ricker model to 4 sites of data.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data
df <- readRDS("data_working/df_207sites.rds")
test2 <- df[c(33,43)]
#df <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_4sites/df_4sites.rds")

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

test_data2 <- lapply(test2, function(x) stan_data_compile(x))
#stan_data_l <- lapply(df, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model

# sets initial values of c and s to help chains converge
init_Ricker <- function(...) {
  list(c = 0.5, s = 0.5)
}

## export results

PM_outputlist_test2 <- lapply(test_data2,
                               function(x) stan("code/teton_moresites/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP.stan",
                                                data = x, chains = 1,iter = 5000,
                                                init = init_Ricker,
                                                control = list(max_treedepth = 12)))

#saveRDS(PM_outputlist_Ricker, "/project/modelscape/users/hlowman/jobresults/teton_4sites/teton_4rivers_output_Ricker_2022_01_22.rds")

# End of script.
