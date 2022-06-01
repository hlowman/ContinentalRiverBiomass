## Fitting models to data
## Step FOUR in Metabolism Modeling Workflow
## June 1, 2022
## Heili Lowman

# I'll be sending a few more jobs to Teton to try and figure out the best
# configuration for the random effects structure.

# This is the script sent today to try the old formulation of parameters,
# but with truncation and predGPP removed. I've also made all parameters
# exponential in form.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data - all in matrix format for new random effects structure
line_lengths <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/line_lengths206.rds")
light_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/light_test10.rds")
gpp_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/gpp_test10.rds")
gpp_sd_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/gpp_sd_test10.rds")
tQ_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/tQ_test10.rds")
e_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_re/e_test10.rds")

rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data for partial run of 10 sites (max of 2374 observations)
stan_data_l <- list(sites = 10, # number of sites 
                    Nobs = 2374, # max number of observations (days)
                    Ndays = line_lengths[c(7:9, 12, 15, 16, 22, 23, 25, 32),], # number of observations per site
                    light = light_mx, # standardized light data
                    GPP = gpp_mx, # standardized GPP estimates
                    GPP_sd = gpp_sd_mx, # standardized GPP standard deviations
                    tQ = tQ_mx, # 10 yr flood standardized discharge
                    new_e = e_mx) # indices to reinitialize biomass estimation

# sets initial values of c and s to help chains converge
init_Ricker <- function(...) {
  list(csite = rep(0.5,length.out = 10), 
       ssite = rep(1.5,length.out = 10),
       rsite = rep(0.3,length.out = 10),
       lsite = rep(-0.05,length.out = 10),
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
PM_outputlist_Ricker_test <- stan("/project/modelscape/users/hlowman/jobscripts/teton_re/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP_re_newpriors_nopredGPP_exp.stan",
                                  data = stan_data_l,
                                  chains = 3,
                                  iter = 5000,
                                  init = init_Ricker,
                                  control = list(max_treedepth = 12))

# export data
saveRDS(PM_outputlist_Ricker_test, "/project/modelscape/users/hlowman/jobresults/teton_re/stan_10rivers_output_Ricker_re_2022_06_01_2.rds")

# End of script.
