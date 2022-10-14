## Resilience of Stream Productivity to Disturbance
## October 13, 2022
## Heili Lowman

# The following set of scripts will walk through the steps necessary to
# prep and send data to Teton as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_raw" and "data_working" folders have been ignored
# using git.ignore, so links to the raw data sets are provided in the step1
# file. If you are accessing the code via GitHub, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running
# the code.

#### Setup ####

# Load necessary packages.
lapply(c("tidyverse", "lubridate", "data.table",
         "rstan","bayesplot","shinystan", "here",
         "GGally", "glmmTMB", "MuMIn", "effects",
         "DHARMa", "lme4", "multcomp", "patchwork",
         "calecopal", "viridis"), require, character.only=T)

# Load necessary datasets.
# Load site-level info (hypoxia and Appling datasets).
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
site <- fread("data_raw/site_data.tsv")

# Load in original list of data fed into the model.
dat_in <- readRDS("data_working/list_182sites_Qmaxnorm_allSL.rds")

# Load in model diagnostics for site-level parameters.
dat_diag <- readRDS("data_working/teton_182rivers_model_diags_101422.rds")

# Load in list containing all iterations of site-level parameters.
dat_out <- readRDS("data_working/teton_182rivers_model_params_all_iterations_101422.rds")

#### Data Prep ####

# Take list containing all input data and make into a df.
dat_in_df <- map_df(dat_in, ~as.data.frame(.x), .id="site_name")

# Take list containing all iterations of parameters and make into a df.
dat_out_df <- map_df(dat_out, ~as.data.frame(.x), .id="site_name")

#### Value filter for rmax ####

# Negative rmax values are not biologically reasonable, so I've 
# removed them.

# First, calculate mean rmax values at all the sites.
dat_out_rmean <- dat_out_df %>%
  group_by(site_name) %>%
  summarize(r_mean = mean(r)) %>%
  ungroup()

# And remove negative values.
dat_out_rmean_pos <- dat_out_rmean %>%
  filter(r_mean > 0) # Removes 7 sites.

#### Rhat filter for rmax ####

# Before proceeding with the first step on my analyses, I will be filtering out 
# sites at which the model did not converge well for the rmax parameter.
# Sites with Rhat > 1.05 will not pass muster.

dat_diag_rfilter1 <- dat_diag %>%
  filter(parameter == "r") %>%
  filter(Rhat < 1.05) # 6 sites drop off

#### normRMSE filter ####

# Next, I will be filtering out sites that do not do a good job of predicting
# GPP using the original data.

# Using code from Joanna's scripts "Predicted_ProductivityModel_Ricker.R"
# and "Biomass2_WSpredictions.R".

# First, need to create the function for predicting GPP.
PM_Ricker <- function(r, lambda, s, c, sig_p, sig_o, df) {
  
  ## Data
  Ndays <- length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$Q_rel # discharge standardized to max value
  new_e <- df$new_e
  
  ## Vectors for model output of P, B, pred_GPP
  P <- numeric(Ndays)
  P[1] <- 1
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*100*(tQ[i] - c)))
  }
  
  B <- numeric(Ndays)
  B[1] <- log(GPP[1]/light[1])
  pred_GPP <- numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])
  
  ## Process Model
  for(j in 2:Ndays){
    # adding in section for my re-initialization functionality
    if (new_e[j]==1) {
      
      B[j] ~ MCMCglmm::rtnorm(1, mean = log(GPP[j]/light[j])) }
    
    else {
      
      B[j] <- MCMCglmm::rtnorm(1, mean = (B[j-1] + r + lambda*exp(B[j-1]))*P[j],
                               sd = sig_p, upper = 5) }
    
  }
  
  for(i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), 
                                    sd = sig_o, lower = 0.01)
  }
  
  return(pred_GPP)
}

# Next, need to write the function with which to perform the simulation.
Ricker_sim_fxn <- function(y, x){
  # identify data
  output <- y # Teton/stan output
  df <- x # original data input
  
  # extracted parameters from STAN output already
  pars <- output
  
  # create empty matrix with days of GPP x length of iterations to receive values
  simmat <- matrix(NA, length(df$GPP), length(unlist(pars$sig_p)))
  rmsemat <- matrix(NA, length(df$GPP), 1)
  
  # simulate pred_GPP holding a parameter set for a given iteration constant
  # and then predict forward for the time period of interest (i.e., length(df$GPP))
  for(i in 1:length(pars$r)){
    simmat[,i] <- PM_Ricker(pars$r[i], pars$lambda[i], pars$s[i], pars$c[i], pars$sig_p[i], pars$sig_p[i], df)
    rmsemat[i] <- sqrt(sum((simmat[,i] - df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat, rmsemat)
  return(l)
  
}

# And finally, apply the function to my data.
Ricker_sim_60sites <- mapply(Ricker_sim_fxn, dat_out, dat_in)

#predGPP <- Ricker_sim_2sites[[1]]
#rmse <- Ricker_sim_2sites[[2]]

# Adding the nRMSE calculation into the function above didn't play nicely with
# the list that existed, so calculating outside instead.
nRMSE_fxn <- function(df, df_orig){
  
  # Calculate the mean RMSE value for each site.
  nRMSE <- mean(df)/(max(df_orig$GPP) - min(df_orig$GPP))
  
}

nRMSE_60sites <- mapply(nRMSE_fxn, Ricker_sim_60sites, dat_in)

nRMSE_60sitesdf <- as.data.frame(nRMSE_60sites) %>%
  rownames_to_column("site_name") %>%
  rename("nRMSE" = "nRMSE_60sites")

# Export both sets of results.
#saveRDS(Ricker_sim_60sites, "data_working/Sim_Ricker_60sites_101322.rds")
saveRDS(nRMSE_60sitesdf, "data_working/nRMSE_60sites_101422.rds")

dat_rmse_rfilter2 <- nRMSE_60sitesdf %>%
  filter(nRMSE < 0.5) #???

#### Figures ####

# Next, append the positive rmax values to the Rhat filter to remove
# appropriate sites.
dat_out_rmean_Rhat <- inner_join(dat_diag_rfilter1, dat_out_rmean_pos)

# Then, remove any sites based on RMSE values.
#dat_out_rmean_Rhat_rmse <- left_join(dat_rmse_rfilter2, dat_out_rmean_Rhat)

# Finally, calculate coefficient of variation in discharge at every site,
# and add to the dataset for plotting purposes.
dat_in_cvq <- dat_in_df %>%
  group_by(site_name) %>%
  summarize(cvQ = (sd(Q, na.rm = TRUE)/mean(Q, na.rm = TRUE))) %>%
  ungroup()

# And, append this to the larger dataset.
dat_out_yas <- left_join(dat_out_rmean_Rhat, dat_in_cvq)

# Distribution of rmax values:
(fig1 <- ggplot(dat_out_yas, aes(x = r_mean)) +
  geom_histogram(bins = 60, alpha = 0.8, 
                 fill = "#A99CD9", color = "#A99CD9") +
  labs(x = "Maximum Growth Rate (r)",
       y = "Count") +
  theme_bw())

# CV of Discharge vs. rmax:
(fig2 <- ggplot(dat_out_yas, aes(x = cvQ, y = r_mean)) +
    geom_point(alpha = 0.8, fill = "808C91") +
    labs(x = "Coefficient of Variation in Discharge (CVq)",
         y = "Maximum Growth Rate (r)") +
    theme_bw())

# Combine figures above.
(fig_r <- fig1 + fig2)

# ggsave(fig_r,
#        filename = "figures/teton_fall22/r_cvQ.jpg",
#        width = 10,
#        height = 5)

# End of script.
