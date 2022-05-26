## Fitting models to data
## Step FOUR in Metabolism Modeling Workflow
## May 25, 2022
## Heili Lowman

# I'll be sending a few jobs to Teton to try and figure out the best
# configuration for the random effects structure.

# This is the second script sent today.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data - all in matrix format for new random effects structure
line_lengths <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/line_lengths206.rds")
light_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/light_mx43.rds")
gpp_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/gpp_mx43.rds")
gpp_sd_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/gpp_sd_mx43.rds")
tQ_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/tQ_mx43.rds")
e_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/e_mx43.rds")

rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data for partial run of 43 sites (max of 2926 observations)
stan_data_l <- list(sites = 43, # number of sites 
                    Nobs = 2926, # max number of observations (days)
                    Ndays = line_lengths[c(7:9, 12, 15, 16, 22, 23, 25, 32, 33, 39, 49, 55, 60, 66, 70:72, 76, 79, 81, 98, 101, 115, 116, 119, 120, 132:134, 140, 142, 157, 163:166, 169, 172, 176, 189, 200),], # number of observations per site
                    light = light_mx, # standardized light data
                    GPP = gpp_mx, # standardized GPP estimates
                    GPP_sd = gpp_sd_mx, # standardized GPP standard deviations
                    tQ = tQ_mx, # 10 yr flood standardized discharge
                    new_e = e_mx) # indices to reinitialize biomass estimation

# sets initial values of c and s to help chains converge
init_Ricker <- function(...) {
  list(csite = rep(0.5,length.out = 43), 
       ssite = rep(1.5,length.out = 43),
       rsite = rep(0.3,length.out = 43),
       lsite = rep(-0.05,length.out = 43),
       c = rep(0.5,length.out = 1), 
       s = rep(1.5,length.out = 1),
       r = rep(0.3,length.out = 1),
       lambda = rep(-0.05,length.out = 1),
       csigma = rep(0.1,length.out = 1), 
       ssigma = rep(0.1,length.out = 1),
       rsigma = rep(0.1,length.out = 1),
       lsigma = rep(0.1,length.out = 1)) # more new values as of MAY 2022
}

## run test model and export results
PM_outputlist_Ricker_test <- stan("/project/modelscape/users/hlowman/jobscripts/teton_re/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP_re_newpriors_notrunc_nopredGPP.stan",
                                  data = stan_data_l,
                                  chains = 3,
                                  iter = 5000,
                                  init = init_Ricker,
                                  control = list(max_treedepth = 12))

# export data
saveRDS(PM_outputlist_Ricker_test, "/project/modelscape/users/hlowman/jobresults/teton_re/stan_43rivers_output_Ricker_re_2022_05_25.rds")

# End of script.
