## 1 river simulation practice script
## Created: September 15, 2021
## Modified: November 30, 2021
## Heili Lowman

# The following code will do some simulation of data to verify
# the fit/approach of the Ricker model.

# Specifically, this code was adapted from an earlier version to
# simulate gappy time series data so that I can test a model
# combination approach.

# Load packages
lapply(c("calecopal", "cowplot", "parallel",
         "devtools", "viridis", "here",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", 
         "ggrepel", "patchwork", "MCMCglmm",
         "data.table"), require, character.only=T)

#### Data Import & Processing ####

# I've chosen to switch to nwis_01608500, South Branch Potomac River (WV),
# because it has an overall good model fit,
# and has baseline light (PAR surface) and discharge (Q) data.

# Load dataset loaded into Teton and filter for desired site.
data_in <- readRDS("data_working/df_34sites.rds")
data_in_1site <- data_in$nwis_01608500

#### Simulation of GPP data ####

# The following is from Joanna's script "Predicted_ProductivityModel_Ricker.R"
# Data simulation function

PM_Ricker <- function(r, lambda, s, c, sig_p, sig_o, df) { # Need to provide model parameters AND original data.
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value
  
  ## Vectors for model output of P, B, pred_GPP
  P <- numeric(Ndays)
  P[1] <- 1
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*100*(tQ[i] - c)))
  }
  
  B<-numeric(Ndays)
  B[1] <- log(GPP[1]/light[1])
  
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] <- MCMCglmm::rtnorm(1, mean = (B[j-1] + r + lambda*exp(B[j-1]))*P[j], sd = sig_p, upper = 5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), sd = sig_o, lower=0.01)
  }
  
  return(pred_GPP)
}

# For nwis_01608500
# Using means in "data_working/teton_34rivers_model_parameters_090821.rds"
# lambda <- -(0.127943684/8.6776636)
gpp_sim <- PM_Ricker(r = 0.1279, lambda = -0.0147, s = 34.2720, c = 0.2348, 
                     sig_p = 0.2520, sig_o = 0.8697, df = data_in_1site)

# Now, to build the dataset that we'll fit the model to
sim_dat <- data_in_1site %>%
  select(date, GPP_sd, light_rel, tQ) %>%
  mutate(GPP = gpp_sim) # NOTE THE RENAMING

#### Fit Ricker model to simulated GPP data ####

# The remaining code below is copied/modified from the Biomass1_FitModels.R script.

## Source data
df <- sim_dat

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data - number of days of data, light, gpp, discharge
stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), 
               light = x$light_rel, 
               GPP = x$GPP,
               GPP_sd = x$GPP_sd, 
               tQ = x$tQ)
  return(data)
}

stan_data_l <- stan_data_compile(df)

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

## PM 2 - Latent Biomass (Ricker)
# With Persistence Term (P)

# sets initial values of c and s to help chain converge
init_Ricker <- function(...) {
  list(c = 0.5, s = 100)
}

## export results
PM_outputlist_Ricker <- stan("code/pooling_practice/Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                                                data = stan_data_l,chains = 3,iter = 5000,
                                                init = init_Ricker,
                                                control = list(max_treedepth = 12))

saveRDS(PM_outputlist_Ricker, "data_working/simulation_1site_output_Ricker_2021_11_30.rds")

#### Re-extraction of model parameters ####

# Extract the parameters resulting from fitting the simulated data to the model.
sim_params <- extract(PM_outputlist_Ricker, c("r","lambda","s","c",
                                                 "B","P","pred_GPP","sig_p","sig_o"))

# And create a dataframe
sim_params_df <- as.data.frame(sim_params) %>%
  # and add "K" to it, calculating for each individual iteration
  mutate(k = (-1*r)/lambda)

# And now to calculate means by site.
sim_means <- sim_params_df %>%
  summarize(r_mean = mean(r),
            k_mean = mean(k),
            lambda_mean = mean(lambda),
            s_mean = mean(s),
            c_mean = mean(c),
            sigp_mean = mean(sig_p),
            sigo_mean = mean(sig_o))

# Parameter     Original Value  Simulated Output
# r             0.1279            0.1472
# lambda        -0.0147          -0.0156
# s             34.2720         339.2987
# c             0.2348            0.2446

# Final thoughts - looks good except for the s value...
# Per discussion with Joanna on 12/2, this is likely due to the x100
# I have in the code above, but wasn't in the STAN model script.

#### Time Sequence Delineation ####

# Borrowing this code from code/teton_moresites/Data_Availability_figures.R
# Using sim_dat dataframe created above
df <- sim_dat

# create placeholder columns for day to day differences and sequences
df <- df %>%
  mutate(diff_time = 0, seq = 1, new_e = 0)

# create function for application to all sites
# so that we can loop over the separate time sequences for a given site
seqFUN <- function(d){
  
  # calculate the difference from one day to the next
  for(i in 2:nrow(d)){
    d$diff_time[i] = difftime(time1 = d$date[i], time2 = d$date[(i-1)],
                              units = "days")
  }
  
  # delineate sequenced time frames based on day to day differences
  # anything less than 14 day gaps is permissible
  for(i in 2:nrow(d)){
    if(d$diff_time[i] < 14){
      d$seq[i] = d$seq[(i-1)]
    } else {
      d$seq[i] = d$seq[(i-1)] + 1
    }
  }
  
  # add column to delineate changes in events
  for(i in 2:nrow(d)){
    d$new_e[i] = d$seq[i] - d$seq[i-1]
  }
  
  return(d)
  
}

# And now map this to the entire site list.
events_dat <- seqFUN(df)

#### Fit New Ricker Model to Simulated GPP Data ####

## Source data
df <- events_dat

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data - number of days of data, light, gpp, discharge
stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), 
               light = x$light_rel, 
               GPP = x$GPP,
               GPP_sd = x$GPP_sd, 
               tQ = x$tQ,
               new_e = x$new_e)
  return(data)
}

# Need to keep this as a list to iterate over each event
stan_data_l <- stan_data_compile(df)

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

## PM 2 - Latent Biomass (Ricker)
# With Persistence Term (P)

# sets initial values of c and s to help chain converge
init_Ricker <- function(...) {
  list(c = 0.5, s = 100)
}

## export results
PM_outputlist_Ricker <- stan("code/pooling_practice/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts.stan",
                             data = stan_data_l,chains = 3,iter = 5000,
                             init = init_Ricker,
                             control = list(max_treedepth = 12))

saveRDS(PM_outputlist_Ricker, "data_working/simulation_1site_output_Ricker_2021_12_02.rds")

#### Re-re-extraction of model parameters ####

# Extract the parameters resulting from fitting the simulated data to the model.
sim_params <- extract(PM_outputlist_Ricker, c("r","lambda","s","c",
                                              "B","P","pred_GPP","sig_p","sig_o"))

# And create a dataframe
sim_params_df <- as.data.frame(sim_params) %>%
  # and add "K" to it, calculating for each individual iteration
  mutate(k = (-1*r)/lambda)

# And now to calculate means by site.
sim_means <- sim_params_df %>%
  summarize(r_mean = mean(r),
            k_mean = mean(k),
            lambda_mean = mean(lambda),
            s_mean = mean(s),
            c_mean = mean(c),
            sigp_mean = mean(sig_p),
            sigo_mean = mean(sig_o))

# Parameter     Original Value  Simulated Output  Simulated Output(without reinit)
# r             0.1279            0.1215            0.1472
# lambda        -0.0147          -0.0143           -0.0156
# s             34.2720         197.0388          339.2987
# c             0.2348            0.2489            0.2446

# Final thoughts - similar results as above!

# End of script.
