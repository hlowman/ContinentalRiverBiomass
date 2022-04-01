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


####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data
stan_data_compile <- function(x){
  data <- list(Ndays = x$Ndays, 
               light = x$light_rel, 
               GPP = x$GPP,
               GPP_sd = x$GPP_sd, 
               tQ = x$tQ,
               new_e = x$new_e) # new column for reinitialization of B[j] values
  
  return(data)
}

stan_data_l <- lapply(df, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model

# sets initial values of c and s to help chains converge
# init_Ricker <- function(...) {
#   list(c = 0.5, s = 0.5) # new values as of jan 2022
# }

## export results

PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("/project/modelscape/users/hlowman/jobscripts/teton_moresites/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_noP.stan",
                                                data = x, chains = 3,iter = 5000,
                                                #init = init_Ricker,
                                                control = list(max_treedepth = 12)))

saveRDS(PM_outputlist_Ricker, "/project/modelscape/users/hlowman/jobresults/teton_moresites/teton_131rivers_output_Ricker_2022_02_07.rds")

# End of script.