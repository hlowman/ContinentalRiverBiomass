## Fitting models to data
## Step FOUR in Metabolism Modeling Workflow
## April 20, 2022
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to fit the Ricker model to 206 sites for my JASM poster presentation.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data - with new index number by site
df <- readRDS("data_working/df_206sites_10yrQnorm.rds")

#### Creating Necessary Matrices ####

## Pull data from the assembled list into a matrix to pass to the stan file

# Pull out only light_rel column from each dataframe within the list
names <- "light_rel"
test_light <- lapply(df, "[", names)

# Count the length of each line
length_df <- function(x){
  data <- length(x$light_rel) # need to go two layers in to calculate length of light_column
  return(data)
}

line_lengths <- lapply(test_light, length_df) # apply function to full list

line_lengths <- t(as.data.frame(line_lengths)) # and transpose into a dataframe

# Pad shorter lines with 0 below
test_light[which(line_lengths != max(line_lengths))] <- 
  lapply(test_light[which(line_lengths != max(line_lengths))], function(x){
     # create list of existing light values
    list1 <- x$light_rel
    list1 <- as.data.frame(list1) %>%
      rename(light_rel = list1) 
     # create list of NAs to be added
    list2 <- rep(0, times = max(line_lengths)-length(x$light_rel))
    list2 <- as.data.frame(list2) %>%
      rename(light_rel = list2)
    
     # join the two together
    rbind(list1, list2)
  })

# Use the list created above to create a matrix of light values
test_light_mx <- matrix(NA, 3208, 206)
test_light_mx <- matrix(unlist(test_light), nrow = 3208, ncol = 206)
# 206 sites
# 3208 days is the longest time series
# Yippee!! This works!!

# Renaming matrix for use in model test below.
light_mx <- test_light_mx

# Repeat for GPP
# Pull out only GPP column from each dataframe within the list
names2 <- "GPP"
subset_gpp <- lapply(df, "[", names2)

# Count the length of each line
length_df2 <- function(x){
  data <- length(x$GPP) # need to go two layers in to calculate length of GPP column
  return(data)
}

line_lengths2 <- lapply(subset_gpp, length_df2) # apply function to full list

line_lengths2 <- t(as.data.frame(line_lengths2)) # and transpose into a dataframe

# Pad shorter lines with 0 below
subset_gpp[which(line_lengths2 != max(line_lengths2))] <- 
  lapply(subset_gpp[which(line_lengths2 != max(line_lengths2))], function(x){
    # create list of existing GPP values
    list1 <- x$GPP
    list1 <- as.data.frame(list1) %>%
      rename(GPP = list1) 
    # create list of NAs to be added
    list2 <- rep(0, times = max(line_lengths2)-length(x$GPP))
    list2 <- as.data.frame(list2) %>%
      rename(GPP = list2)
    # join the two together
    rbind(list1, list2)
  })

# Use the list created above to create a matrix of GPP values
gpp_mx <- matrix(NA, 3208, 206)
gpp_mx <- matrix(unlist(subset_gpp), nrow = 3208, ncol = 206)
# WOOHOO

# And for sd of GPP.
# Pull out only GPP_sd column from each dataframe within the list
names3 <- "GPP_sd"
subset_gppsd <- lapply(df, "[", names3)

# Count the length of each line
length_df3 <- function(x){
  data <- length(x$GPP_sd) # need to go two layers in to calculate length of GPP_sd column
  return(data)
}

line_lengths3 <- lapply(subset_gppsd, length_df3) # apply function to full list

line_lengths3 <- t(as.data.frame(line_lengths3)) # and transpose into a dataframe

# Pad shorter lines with 0 below
subset_gppsd[which(line_lengths3 != max(line_lengths3))] <- 
  lapply(subset_gppsd[which(line_lengths3 != max(line_lengths3))], function(x){
    # create list of existing GPP sd values
    list1 <- x$GPP_sd
    list1 <- as.data.frame(list1) %>%
      rename(GPP_sd = list1) 
    # create list of NAs to be added
    list2 <- rep(0, times = max(line_lengths3)-length(x$GPP_sd))
    list2 <- as.data.frame(list2) %>%
      rename(GPP_sd = list2)
    # join the two together
    rbind(list1, list2)
  })

# Use the list created above to create a matrix of GPP values
gpp_sd_mx <- matrix(NA, 3208, 206)
gpp_sd_mx <- matrix(unlist(subset_gppsd), nrow = 3208, ncol = 206)
# WOOT

# And for discharge.
# Pull out only tQ column from each dataframe within the list
names4 <- "tQ"
subset_tQ <- lapply(df, "[", names4)

# Count the length of each line
length_df4 <- function(x){
  data <- length(x$tQ) # need to go two layers in to calculate length of tQ column
  return(data)
}

line_lengths4 <- lapply(subset_tQ, length_df4) # apply function to full list

line_lengths4 <- t(as.data.frame(line_lengths4)) # and transpose into a dataframe

# Pad shorter lines with 0 below
subset_tQ[which(line_lengths4 != max(line_lengths4))] <- 
  lapply(subset_tQ[which(line_lengths4 != max(line_lengths4))], function(x){
    # create list of existing tQ values
    list1 <- x$tQ
    list1 <- as.data.frame(list1) %>%
      rename(tQ = list1) 
    # create list of NAs to be added
    list2 <- rep(0, times = max(line_lengths4)-length(x$tQ))
    list2 <- as.data.frame(list2) %>%
      rename(tQ = list2)
    # join the two together
    rbind(list1, list2)
  })

# Use the list created above to create a matrix of tQ values
tQ_mx <- matrix(NA, 3208, 206)
tQ_mx <- matrix(unlist(subset_tQ), nrow = 3208, ncol = 206)

# and finally the reinitialization index values.
# Pull out only new_e column from each dataframe within the list
names5 <- "new_e"
subset_e <- lapply(df, "[", names5)

# Count the length of each line
length_df5 <- function(x){
  data <- length(x$new_e) # need to go two layers in to calculate length of new_e column
  return(data)
}

line_lengths5 <- lapply(subset_e, length_df5) # apply function to full list

line_lengths5 <- t(as.data.frame(line_lengths5)) # and transpose into a dataframe

# Pad shorter lines with 0 below
subset_e[which(line_lengths5 != max(line_lengths5))] <- 
  lapply(subset_e[which(line_lengths5 != max(line_lengths5))], function(x){
    # create list of existing e values
    list1 <- x$new_e
    list1 <- as.data.frame(list1) %>%
      rename(new_e = list1) 
    # create list of NAs to be added
    list2 <- rep(0, times = max(line_lengths5)-length(x$new_e))
    list2 <- as.data.frame(list2) %>%
      rename(new_e = list2)
    # join the two together
    rbind(list1, list2)
  })

# Use the list created above to create a matrix of tQ values
e_mx <- matrix(NA, 3208, 206)
e_mx <- matrix(unlist(subset_e), nrow = 3208, ncol = 206)

# Export these matrices for use in Teton runs.
# Note - re-exported these on May 4, 2022 since I re-normalized
# discharge values using 10 year flood in cms, not cfs.
saveRDS(line_lengths, "data_working/line_lengths206.rds")
saveRDS(light_mx, "data_working/light_mx206.rds")
saveRDS(gpp_mx, "data_working/gpp_mx206.rds")
saveRDS(gpp_sd_mx, "data_working/gpp_sd_mx206.rds")
saveRDS(tQ_mx, "data_working/tQ_mx206.rds")
saveRDS(e_mx, "data_working/e_mx206.rds")

# Subset to two sites of exact same length for test run prior to sending full job to Teton.
light_mx2 <- light_mx[1:192,c(3:4)]
gpp_mx2 <- gpp_mx[1:192,c(3:4)]
gpp_sd_mx2 <- gpp_sd_mx[1:192,c(3:4)]
tQ_mx2 <- tQ_mx[1:192,c(3:4)]
e_mx2 <- e_mx[1:192,c(3:4)]

# Subset to ten sites of a similar length for test run prior to sending full job to Teton.
light_mx10 <- light_mx[1:192,c(1:5, 31, 40, 78, 94, 129)]
gpp_mx10 <- gpp_mx[1:192,c(1:5, 31, 40, 78, 94, 129)]
gpp_sd_mx10 <- gpp_sd_mx[1:192,c(1:5, 31, 40, 78, 94, 129)]
tQ_mx10 <- tQ_mx[1:192,c(1:5, 31, 40, 78, 94, 129)]
e_mx10 <- e_mx[1:192,c(1:5, 31, 40, 78, 94, 129)]

# Ok, so server isn't loading R for some reason, so proceeding on desktop.

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data for 2 site test run (max 192 observations)
stan_data_l2 <- list(sites = 2, # number of sites 
                      Nobs = 192, # max number of observations (days)
                      Ndays = line_lengths[c(3:4),], # number of observations per site
                      light = light_mx2, # standardized light data
                      GPP = gpp_mx2, # standardized GPP estimates
                      GPP_sd = gpp_sd_mx2, # standardized GPP standard deviations
                      tQ = tQ_mx2, # 10 yr flood standardized discharge
                      new_e = e_mx2) # indices denoting when to reinitialize biomass estimation

## compile data for 10 site test run (max 192 observations)
stan_data_l10 <- list(sites = 10, # number of sites 
                    Nobs = 192, # max number of observations (days)
                    Ndays = line_lengths[c(1:5, 31, 40, 78, 94, 129),], # number of observations per site
                    light = light_mx10, # standardized light data
                    GPP = gpp_mx10, # standardized GPP estimates
                    GPP_sd = gpp_sd_mx10, # standardized GPP standard deviations
                    tQ = tQ_mx10, # 10 yr flood standardized discharge
                    new_e = e_mx10) # indices denoting when to reinitialize biomass estimation

saveRDS(stan_data_l10, "data_working/stan_10rivers_input_list.rds")

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model

# sets initial values to help chains converge
# this was particularly problematic for 10+ sites
init_Ricker <- function(...) {
  list(csite = rep(0.5,length.out = 10), 
       ssite = rep(0.5,length.out = 10),
       rsite = rep(0.5,length.out = 10),
       lsite = rep(0.5,length.out = 10)) # new values as of apr 2022
}

# PM_outputlist_Ricker2 <- stan("code/teton_jasm/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP_re.stan",
#                                data = stan_data_l2, chains = 1,iter = 5000,
#                                init = init_Ricker,
#                                control = list(max_treedepth = 12)) # works

# The following takes ~3 hours to run on desktop.
PM_outputlist_Ricker10 <- stan("code/teton_jasm/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP_re.stan",
                            data = stan_data_l10, chains = 3,iter = 5000,
                            init = init_Ricker,
                            control = list(max_treedepth = 12)) # RUNS!! OMGEEEE!!!

# At first, it was breaking but the message output during the
# initialization step is the same as for some runs that do work.
# Chain 1: Rejecting initial value:
# Chain 1:   Log probability evaluates to log(0), i.e. negative infinity.
# Chain 1:   Stan can't start sampling from this initial value.

# However, this was fixed by giving the sampler an initialization value of 0.5
# for all four parameters (csite, ssite, rsite, lsite).

launch_shinystan(PM_outputlist_Ricker10)
# Well, it looks terrible, but it ran.

## export results
saveRDS(PM_outputlist_Ricker10, "data_working/stan_10rivers_output_Ricker_re_2022_04_22.rds")

#### Additional runs performed on the server ####

# In response to the terrible number of divergent transitions, I'm going to run
# a 10 site model on the server to speed things up but try and change the inits
# per Joanna's suggestion.

# First, make sure I'm using the matrices from above.
## Source data - all in matrix format for new random effects structure
line_lengths <- readRDS("data_working/line_lengths206.rds")
light_mx <- readRDS("data_working/light_mx206.rds")
gpp_mx <- readRDS("data_working/gpp_mx206.rds")
gpp_sd_mx <- readRDS("data_working/gpp_sd_mx206.rds")
tQ_mx <- readRDS("data_working/tQ_mx206.rds")
e_mx <- readRDS("data_working/e_mx206.rds")

# And trim down to first 10 good sites
light_test10_again <- light_mx[1:2374,c(7:9, 12, 15, 16, 22, 23, 25, 32)]
gpp_test10_again <- gpp_mx[1:2374,c(7:9, 12, 15, 16, 22, 23, 25, 32)]
gpp_sd_test10_again <- gpp_sd_mx[1:2374,c(7:9, 12, 15, 16, 22, 23, 25, 32)]
tQ_test10_again <- tQ_mx[1:2374,c(7:9, 12, 15, 16, 22, 23, 25, 32)]
e_test10_again <- e_mx[1:2374,c(7:9, 12, 15, 16, 22, 23, 25, 32)]

rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data for partial run of 10 sites (max of 2374 observations)
stan_data_l_test10_again <- list(sites = 10, # number of sites 
                    Nobs = 2374, # max number of observations (days)
                    Ndays = line_lengths[c(7:9, 12, 15, 16, 22, 23, 25, 32),], # number of observations per site
                    light = light_test10_again, # standardized light data
                    GPP = gpp_test10_again, # standardized GPP estimates
                    GPP_sd = gpp_sd_test10_again, # standardized GPP standard deviations
                    tQ = tQ_test10_again, # 10 yr flood standardized discharge
                    new_e = e_test10_again) # indices denoting when to reinitialize biomass estimation

# sets initial values of c and s to help chains converge
init_Ricker <- function(...) {
  list(csite = rep(0.25,length.out = 10), 
       ssite = rep(1.5,length.out = 10),
       rsite = rep(0.3,length.out = 10),
       lsite = rep(-0.05,length.out = 10)) # new values as of MAY 2022
}

## run test model and export results
# started at 10:33AM
PM_outputlist_Ricker_test10_again <- stan("code/teton_jasm/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP_re_newpriors.stan",
                                      data = stan_data_l_test10_again,
                                      chains = 3,
                                      iter = 5000,
                                      init = init_Ricker,
                                      control = list(max_treedepth = 12))

# Error message displayed:
# Warning messages:
#   1: In validityMethod(object) :
#   The following variables have undefined values:  pred_GPP[1197,1],The following variables have undefined values:  pred_GPP[1198,1],The following variables have undefined values:  pred_GPP[1199,1],The following variables have undefined values:  pred_GPP[1200,1],The following variables have undefined values:  pred_GPP[1201,1],The following variables have undefined values:  pred_GPP[1202,1],The following variables have undefined values:  pred_GPP[1203,1],The following variables have undefined values:  pred_GPP[1204,1],The following variables have undefined values:  pred_GPP[1205,1],The following variables have undefined values:  pred_GPP[1206,1],The following variables have undefined values:  pred_GPP[1207,1],The following variables have undefined values:  pred_GPP[1208,1],The following variables have undefined values:  pred_GPP[1209,1],The following variables have undefined values:  pred_GPP[1210,1],The following variables have undefined values:  pred_GPP[1211,1],The following variables h [... truncated]
# 2: There were 7500 divergent transitions after warmup. See
# https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 3: Examine the pairs() plot to diagnose sampling problems
# 
# 4: The largest R-hat is NA, indicating chains have not mixed.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#r-hat 
# 5: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#bulk-ess 
# 6: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#tail-ess 

launch_shinystan(PM_outputlist_Ricker_test10_again)
# Well, ...

saveRDS(PM_outputlist_Ricker_test10_again, "data_working/stan_10rivers_output_Ricker_re_2022_05_04.rds")

#### Models run after chatting with Joanna and Bob 5/10/22 ####

# Attempt #1: added to initialization values and removed erroneous rsigma and lsigma from 
# STAN script

rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data for partial run of 10 sites (max of 2374 observations)
stan_data_l_test10_again_again <- list(sites = 10, # number of sites 
                                 Nobs = 2374, # max number of observations (days)
                                 Ndays = line_lengths[c(7:9, 12, 15, 16, 22, 23, 25, 32),], # number of observations per site
                                 light = light_test10_again, # standardized light data
                                 GPP = gpp_test10_again, # standardized GPP estimates
                                 GPP_sd = gpp_sd_test10_again, # standardized GPP standard deviations
                                 tQ = tQ_test10_again, # 10 yr flood standardized discharge
                                 new_e = e_test10_again) # indices denoting when to reinitialize biomass estimation

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
# started at 5:13PM
PM_outputlist_Ricker_test10_again_again <- stan("code/teton_jasm/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP_re_newpriors.stan",
                                          data = stan_data_l_test10_again_again,
                                          chains = 3,
                                          iter = 5000,
                                          init = init_Ricker,
                                          control = list(max_treedepth = 12))

launch_shinystan(PM_outputlist_Ricker_test10_again_again)
# Still looks terrible.

saveRDS(PM_outputlist_Ricker_test10_again_again, "data_working/stan_10rivers_output_Ricker_re_2022_05_10.rds")

# End of script.
