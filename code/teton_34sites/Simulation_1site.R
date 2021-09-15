## 1 river simulation practice script
## Created: September 15, 2021
## Heili Lowman

# The following code will do some simulation of data to verify
# the fit/approach of the Ricker model.

# Load packages
lapply(c("calecopal", "cowplot", "parallel",
         "devtools", "viridis", "here",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", 
         "ggrepel", "patchwork", "MCMCglmm"), require, character.only=T)

#### Data Import & Processing ####

# I've chosen to use nwis_03067510, because it has an overall good model fit,
# and has baseline light (PAR surface) and discharge (Q) data.

# Load dataset loaded into Teton and filter for desired site.
data_in <- readRDS("data_working/df_34sites.rds")
data_in_1site <- data_in$nwis_03067510

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

# For nwis_03067510
# Using means in "data_working/teton_34rivers_model_parameters_090821.rds"
gpp_sim <- PM_Ricker(r = 0.02409658, lambda = -0.0206386, s = 109.506, c = 1.167247, sig_p = 0.2728436, sig_o = 0.1573847, df = data_in_1site)

# Now, to build the dataset that we'll fit the model to
sim_dat <- data_in_1site %>%
  select(date, GPP_sd, light_rel, tQ) %>%
  mutate(GPP = gpp_sim)

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

## compile data
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
PM_outputlist_Ricker <- stan("code/teton_34sites/Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                                                data = stan_data_l,chains = 3,iter = 5000,
                                                init = init_Ricker,
                                                control = list(max_treedepth = 12))

saveRDS(PM_outputlist_Ricker, "data_working/simulation_1site_output_Ricker_2021_09_15.rds")

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

# Parameter     Original Value     Simulated Output
# r             0.02409658         0.03597256
# lambda        -0.0206386         -0.04077765
# s             109.506            164.0564
# c             1.167247           1.551805

# End of script.
