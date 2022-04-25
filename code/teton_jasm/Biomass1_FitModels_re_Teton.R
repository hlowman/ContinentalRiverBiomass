## Fitting models to data
## Step FOUR in Metabolism Modeling Workflow
## April 22, 2022
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to fit the Ricker model to 206 sites for my JASM poster presentation.

# This script is actually the one run on Teton.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data - all in matrix format for new random effects structure
line_lengths <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/line_lengths206.rds")
light_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/light_mx206.rds")
gpp_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/gpp_mx206.rds")
gpp_sd_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/gpp_sd_mx206.rds")
tQ_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/tQ_mx206.rds")
e_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/e_mx206.rds")

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data for full run of 206 sites (max of 3208 observations)
stan_data_l <- list(sites = 206, # number of sites 
                    Nobs = 3208, # max number of observations (days)
                    Ndays = line_lengths[c(1:206),], # number of observations per site
                    light = light_mx, # standardized light data
                    GPP = gpp_mx, # standardized GPP estimates
                    GPP_sd = gpp_sd_mx, # standardized GPP standard deviations
                    tQ = tQ_mx, # 10 yr flood standardized discharge
                    new_e = e_mx) # indices denoting when to reinitialize biomass estimation

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model

# sets initial values of c and s to help chains converge
init_Ricker <- function(...) {
  list(csite = rep(0.5,length.out = 206), 
       ssite = rep(0.5,length.out = 206),
       rsite = rep(0.5,length.out = 206),
       lsite = rep(0.5,length.out = 206)) # new values as of apr 2022
}

## export results

PM_outputlist_Ricker <- stan("/project/modelscape/users/hlowman/jobscripts/teton_jasm/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP_re.stan",
                             data = stan_data_l, 
                             chains = 3,
                             iter = 5000,
                             init = init_Ricker,
                             control = list(max_treedepth = 12))

saveRDS(PM_outputlist_Ricker, "/project/modelscape/users/hlowman/jobresults/teton_jasm/teton_206rivers_output_Ricker_2022_04_24.rds")

# End of script.
