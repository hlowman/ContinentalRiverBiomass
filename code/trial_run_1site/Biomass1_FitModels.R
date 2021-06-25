## Fitting models to data
## JR Blaszczak
## June 23, 2021
## Heili Lowman

# I'll be modifying some of Joanna's code from the RiverBiomass repository
# to practice fitting the Ricker model to 1 year of data at a "good" and
# "bad" site.

# I've commented out those steps that I feel, for the time being, I don't
# need to perform, and I've changed the appropriate filepaths to match my
# repository structure.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","Metrics","MCMCglmm","tictoc",
         "here"), require, character.only=T)

## Source data - sources the file itself
source("code/trial_run_1site/DataSource_6rivers_StreamLight.R")

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)#parallel::detectCores())

## compile data
stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ)
  return(data)
}

stan_data_l <- lapply(df, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates - initial tests
#########################################

## Initial tests
#AR
# test_ar <- stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
#              data=stan_data_l$nwis_01608500,
#              chains=3,iter=5000, control=list(max_treedepth=12))
# launch_shinystan(test_ar)

#Ricker

#Ricker - time varying r
init_Ricker <- function(...) {
  list(c = 0.5, s = 100)
}

## initial Ricker model run
test_ricker <- stan("code/trial_run_1site/Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                    data=stan_data_l$nwis_05406457,
                    init = init_Ricker,
                    chains=3,iter=5000, control=list(max_treedepth=12))
launch_shinystan(test_ricker)

#Ricker - time varying r
init_Ricker <- function(...) {
  list(c = 0.5, s = 100)
}
test_ricker_tvr <- stan("Stan_ProductivityModel2_Ricker_fixedinit_r.stan",
             data=stan_data_l$nwis_05406457,
             init = init_Ricker,
             chains=3,iter=5000, control=list(max_treedepth=12))

## examine model - it didn't converge :(
launch_shinystan(test_ricker_tvr)


#Gompertz
# init_Gompertz <- function(...) {
#   list(c = 0.5, s = 100)
# }
# test_Gompertz <- stan("Stan_ProductivityModel3_Gompertz_fixedinit_obserr.stan",
#                     data=stan_data_l$nwis_01608500,
#                     init = init_Gompertz,
#                     chains=3,iter=5000, control=list(max_treedepth=12))
# launch_shinystan(test_Gompertz)


#########################################
## Run Stan to get parameter estimates - all sites
#########################################

## PM 1 - Standard time series
# PM_outputlist_AR <- lapply(stan_data_l,
#                            function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
#                                                    data=x,chains=3,iter=5000, control=list(max_treedepth=12)))
# PM_AR_elapsedtime <- lapply(PM_outputlist_AR, function(x) return(get_elapsed_time(x)))
# saveRDS(PM_outputlist_AR, "./rds files/stan_6riv_output_AR_2021_06_01.rds")
# saveRDS(PM_AR_elapsedtime, "./rds files/stan_6riv_AR_time_2021_06_01.rds")

## PM 2 - Latent Biomass (Ricker)
init_Ricker <- function(...) {
  list(c = 0.5, s = 100)
}

## export results
PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("code/trial_run_1site/Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                                                data=x,chains=3,iter=10000,init = init_Ricker,
                                                control=list(max_treedepth=12)))

# Ran with 5000 iterations
# Warning messages:
#   1: There were 3451 divergent transitions after warmup. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 2: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
# http://mc-stan.org/misc/warnings.html#bfmi-low 
# 3: Examine the pairs() plot to diagnose sampling problems
# 
# 4: The largest R-hat is 2.55, indicating chains have not mixed.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#r-hat 
# 5: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#bulk-ess 
# 6: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#tail-ess

# Ran again with 10000 iterations
# Warning messages:
#   1: There were 1363 divergent transitions after warmup. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 2: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
# http://mc-stan.org/misc/warnings.html#bfmi-low 
# 3: Examine the pairs() plot to diagnose sampling problems

PM_Ricker_elapsedtime <- lapply(PM_outputlist_Ricker, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_Ricker, "data_working/stan_1riv_output_Ricker_10000it_2021_06_23.rds")
saveRDS(PM_Ricker_elapsedtime, "data_working/stan_1riv_Ricker_time_10000it_2021_06_23.rds")

## PM 2 - Latent Biomass (Ricker) - time varying r
# init_Ricker <- function(...) {
#   list(c = 0.5, s = 100)
# }
# 
# PM_outputlist_Ricker_tvr <- lapply(stan_data_l,
#                                function(x) stan("Stan_ProductivityModel2_Ricker_fixedinit_r.stan",
#                                                 data=x,chains=3,iter=5000,init = init_Ricker,
#                                                 control=list(max_treedepth=12)))
# PM_Ricker_elapsedtime_tvr <- lapply(PM_outputlist_Ricker_tvr, function(x) return(get_elapsed_time(x)))
# saveRDS(PM_outputlist_Ricker_tvr, "./rds files/stan_6riv_output_Ricker_tvr_2021_06_01.rds")
# saveRDS(PM_Ricker_elapsedtime_tvr, "./rds files/stan_6riv_Ricker_tvr_time_2021_06_01.rds")

## PM 3 - Latent Biomass (Gompertz)
# init_Gompertz <- function(...) {
#   list(c = 0.5, s = 200)
# }
# 
# PM_outputlist_Gompertz <- lapply(stan_data_l,
#                                  function(x) stan("Stan_ProductivityModel3_Gompertz_fixedinit_obserr.stan",
#                                                   data=x,chains=3,iter=5000, 
#                                                   control=list(max_treedepth=12)))
# PM_Gompertz_elapsedtime <- lapply(PM_outputlist_Gompertz, function(x) return(get_elapsed_time(x)))
# 
# saveRDS(PM_outputlist_Gompertz, "./rds files/stan_6riv_output_Gompertz_2021_06_01.rds")
# saveRDS(PM_Gompertz_elapsedtime, "./rds files/stan_6riv_Gompertz_time_2021_06_01.rds")

## View
# launch_shinystan(PM_outputlist_AR$nwis_08447300)
PM_outputlist_Ricker <- readRDS("data_working/stan_1riv_output_Ricker_10000it_2021_06_23.rds")
launch_shinystan(PM_outputlist_Ricker$nwis_05406457)
# launch_shinystan(PM_outputlist_Ricker_tvr$nwis_11044000)

# Extract parameters
mparams<-extract(PM_outputlist_Ricker$nwis_05406457, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o"))

# End of script.
