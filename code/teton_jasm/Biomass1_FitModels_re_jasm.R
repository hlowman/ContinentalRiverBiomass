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

# Export these matrices for use in the actual Teton run.
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
# a 10 site model on the server to speed things up but try and change the allowed
# variance on parameter sigma values per Ana's suggestion.

# First, load in 43 site dataset we know performed poorly.
## Source data - all in matrix format for new random effects structure
line_lengths <- readRDS("data_working/line_lengths206.rds")
light_mx <- readRDS("data_working/light_mx43.rds")
gpp_mx <- readRDS("data_working/gpp_mx43.rds")
gpp_sd_mx <- readRDS("data_working/gpp_sd_mx43.rds")
tQ_mx <- readRDS("data_working/tQ_mx43.rds")
e_mx <- readRDS("data_working/e_mx43.rds")

# And trim down to first 5 sites
light_test5 <- light_mx[1:2374,c(1:5)]
gpp_test5 <- gpp_mx[1:2374,c(1:5)]
gpp_sd_test5 <- gpp_sd_mx[1:2374,c(1:5)]
tQ_test5 <- tQ_mx[1:2374,c(1:5)]
e_test5 <- e_mx[1:2374,c(1:5)]

rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data for partial run of 5 sites (max of 2374 observations)
stan_data_l_test5 <- list(sites = 5, # number of sites 
                    Nobs = 2374, # max number of observations (days)
                    Ndays = line_lengths[c(7:9, 12, 15),], # number of observations per site
                    light = light_test5, # standardized light data
                    GPP = gpp_test5, # standardized GPP estimates
                    GPP_sd = gpp_sd_test5, # standardized GPP standard deviations
                    tQ = tQ_test5, # 10 yr flood standardized discharge
                    new_e = e_test5) # indices denoting when to reinitialize biomass estimation

# sets initial values of c and s to help chains converge
init_Ricker <- function(...) {
  list(csite = rep(0.5,length.out = 5), 
       ssite = rep(0.5,length.out = 5),
       rsite = rep(0.5,length.out = 5),
       lsite = rep(0.5,length.out = 5)) # new values as of apr 2022
}

## export results
# cauchy(0,10) - ten sites
# took too long on the server so i cancelled it
# PM_outputlist_Ricker_test10_1 <- stan("code/teton_jasm/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP_re_wide.stan",
#                              data = stan_data_l_test10, 
#                              chains = 3,
#                              iter = 5000,
#                              init = init_Ricker,
#                              control = list(max_treedepth = 12))

# cauchy(0, 100) - five sites
# also stopped this - the second chain taking FOR EV ER to converge.
# PM_outputlist_Ricker_test5_1 <- stan("code/teton_jasm/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP_re_wide.stan",
#                                       data = stan_data_l_test5, 
#                                       chains = 3,
#                                       iter = 5000,
#                                       init = init_Ricker,
#                                       control = list(max_treedepth = 12))

# End of script.
