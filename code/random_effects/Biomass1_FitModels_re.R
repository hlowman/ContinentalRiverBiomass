## Fitting models to data
## Step FOUR in Metabolism Modeling Workflow
## March 21, 2022
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to fit the Ricker model to a few sites of data to get the random effect working.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data - with new index number by site
df <- readRDS("data_working/df_207sites_indexed.rds")

#### Creating Necessary Matrices ####

# Testing how to pull data from the assembled list into a matrix to pass to the stan file
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

# Pad shorter lines with NA below
test_light[which(line_lengths != max(line_lengths))] <- 
  lapply(test_light[which(line_lengths != max(line_lengths))], function(x){
     # create list of existing light values
    list1 <- x$light_rel
    list1 <- as.data.frame(list1) %>%
      rename(light_rel = list1) 
     # create list of NAs to be added
    list2 <- rep(NA, times = max(line_lengths)-length(x$light_rel))
    list2 <- as.data.frame(list2) %>%
      rename(light_rel = list2)
    
     # join the two together
    rbind(list1, list2)
  })

# Use the list created above to create a matrix of light values
test_light_mx <- matrix(NA, 3208, 207)
test_light_mx <- matrix(unlist(test_light), nrow = 3208, ncol = 207)
# 207 sites
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

# Pad shorter lines with NA below
subset_gpp[which(line_lengths2 != max(line_lengths2))] <- 
  lapply(subset_gpp[which(line_lengths2 != max(line_lengths2))], function(x){
    # create list of existing GPP values
    list1 <- x$GPP
    list1 <- as.data.frame(list1) %>%
      rename(GPP = list1) 
    # create list of NAs to be added
    list2 <- rep(NA, times = max(line_lengths2)-length(x$GPP))
    list2 <- as.data.frame(list2) %>%
      rename(GPP = list2)
    # join the two together
    rbind(list1, list2)
  })

# Use the list created above to create a matrix of GPP values
gpp_mx <- matrix(NA, 3208, 207)
gpp_mx <- matrix(unlist(subset_gpp), nrow = 3208, ncol = 207)
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

# Pad shorter lines with NA below
subset_gppsd[which(line_lengths3 != max(line_lengths3))] <- 
  lapply(subset_gppsd[which(line_lengths3 != max(line_lengths3))], function(x){
    # create list of existing GPP sd values
    list1 <- x$GPP_sd
    list1 <- as.data.frame(list1) %>%
      rename(GPP_sd = list1) 
    # create list of NAs to be added
    list2 <- rep(NA, times = max(line_lengths3)-length(x$GPP_sd))
    list2 <- as.data.frame(list2) %>%
      rename(GPP_sd = list2)
    # join the two together
    rbind(list1, list2)
  })

# Use the list created above to create a matrix of GPP values
gpp_sd_mx <- matrix(NA, 3208, 207)
gpp_sd_mx <- matrix(unlist(subset_gppsd), nrow = 3208, ncol = 207)
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

# Pad shorter lines with NA below
subset_tQ[which(line_lengths4 != max(line_lengths4))] <- 
  lapply(subset_tQ[which(line_lengths4 != max(line_lengths4))], function(x){
    # create list of existing tQ values
    list1 <- x$tQ
    list1 <- as.data.frame(list1) %>%
      rename(tQ = list1) 
    # create list of NAs to be added
    list2 <- rep(NA, times = max(line_lengths4)-length(x$tQ))
    list2 <- as.data.frame(list2) %>%
      rename(tQ = list2)
    # join the two together
    rbind(list1, list2)
  })

# Use the list created above to create a matrix of tQ values
tQ_mx <- matrix(NA, 3208, 207)
tQ_mx <- matrix(unlist(subset_tQ), nrow = 3208, ncol = 207)

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

# Pad shorter lines with NA below
subset_e[which(line_lengths5 != max(line_lengths5))] <- 
  lapply(subset_e[which(line_lengths5 != max(line_lengths5))], function(x){
    # create list of existing e values
    list1 <- x$new_e
    list1 <- as.data.frame(list1) %>%
      rename(new_e = list1) 
    # create list of NAs to be added
    list2 <- rep(NA, times = max(line_lengths5)-length(x$new_e))
    list2 <- as.data.frame(list2) %>%
      rename(new_e = list2)
    # join the two together
    rbind(list1, list2)
  })

# Use the list created above to create a matrix of tQ values
e_mx <- matrix(NA, 3208, 207)
e_mx <- matrix(unlist(subset_e), nrow = 3208, ncol = 207)

# Subset to first two sites for test run.
light_mx2 <- light_mx[1:101,1:2]
gpp_mx2 <- gpp_mx[1:101,1:2]
gpp_sd_mx2 <- gpp_sd_mx[1:101,1:2]
tQ_mx2 <- tQ_mx[1:101,1:2]
e_mx2 <- e_mx[1:101,1:2]

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data
stan_data_l <- list(sites = 207, # number of sites 
               Nobs = 3208, # max number of observations (days)
               Ndays = line_lengths, # number of observations per site
               light = light_mx, # standardized light data
               GPP = gpp_mx, # standardized GPP estimates
               GPP_sd = gpp_sd_mx, # standardized GPP standard deviations
               tQ = tQ_mx, # standardized discharge
               new_e = e_mx) # indices denoting when to reinitialize biomass estimation

## compile data for 2 site test run (max 101 observations)
stan_data_l2 <- list(sites = 2, # number of sites 
                    Nobs = 101, # max number of observations (days)
                    Ndays = line_lengths[1:2,], # number of observations per site
                    light = light_mx2, # standardized light data
                    GPP = gpp_mx2, # standardized GPP estimates
                    GPP_sd = gpp_sd_mx2, # standardized GPP standard deviations
                    tQ = tQ_mx2, # standardized discharge
                    new_e = e_mx2) # indices denoting when to reinitialize biomass estimation


#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model

# sets initial values of c and s to help chains converge
init_Ricker <- function(...) {
  list(c = 0.5, s = 0.5) # new values as of jan 2022
}

## export results

PM_outputlist_Ricker2 <- stan("code/random_effects/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP_re.stan",
                            data = stan_data_l2, chains = 3,iter = 5000,
                            #init = init_Ricker, 
                            control = list(max_treedepth = 12))

launch_shinystan(PM_outputlist_Ricker2)

saveRDS(PM_outputlist_Ricker2, "data_working/stan_2rivers_output_Ricker_re_2022_04_05.rds")

# End of script.