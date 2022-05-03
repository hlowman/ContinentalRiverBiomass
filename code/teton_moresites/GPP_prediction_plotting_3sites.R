## Plotting GPP vs. pred_GPP
## April 30, 2022
## Heili Lowman

# The following script will pull in the results of the previous modeling
# effort to generate plots of modeled vs. predicted GPP for the JASM poster.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data - 
data_in <- readRDS("data_working/df_207sites.rds")
## Teton output by site (because it's too big otherwise) -
data_out_mo <- readRDS("data_working/teton_mo_model_all_params_all_iterations_020622.rds")
data_out_tx <- readRDS("data_working/teton_tx_model_all_params_all_iterations_020622.rds")
data_out_pa <- readRDS("data_working/teton_pa_model_all_params_all_iterations_020622.rds")
## Site info -
data_info <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")
## Tidied parameter data -
data_params <- readRDS("data_working/teton_bothmodels_parameters_means_021322.rds") %>%
  filter(model == "with P") %>%
  filter(r_mean > 0)
## Divergences data - 
data_divs <- readRDS("data_working/teton_207rivers_model_divergences_bothmodels_021322.rds")

# what are roughly the 2.5, 50, and 97.5%iles of rmax values?
low <- quantile(data_params$r_mean, probs = 0.025) # 0.02
med <- quantile(data_params$r_mean, probs = 0.5) # 0.12
high <- quantile(data_params$r_mean, probs = .975) # 0.62

# identify sites where divergences are low and rmax aligns close to the above quantiles
data_divs10 <- data_divs %>%
  filter(model == "with P") %>%
  filter(divergences < 10)

data_best <- left_join(data_divs10, data_params)

# 2.5% - nwis_07061270 - East Fork Black River near Lesterville, MO
# 50% - nwis_07332622 - Bois D'Arc Ck at FM 409 nr Honey Grove, TX
# 97.5% - nwis_03025500 - Allegheny River at Franklin, PA

#### New Prediction of GPP Workflow ####

# Previously, I was just pulling model fit. The below workflow should now
# re-predict GPP values for plotting purposes.

# Using code from Joanna's scripts "Predicted_ProductivityModel_Ricker.R"
# and "Biomass2_WSpredictions.R".

# First, need to create the function for predicting GPP.
PM_Ricker <- function(r, lambda, s, c, sig_p, sig_o, df) {
  
  ## Data
  Ndays <- length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value
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
  #pars <- extract(output, c("r", "lambda", "s", "c", "B", "P", "pred_GPP", "sig_p", "sig_o"))
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
Ricker_sim_mo <- Ricker_sim_fxn(data_mo_site, data_in$nwis_07061270)
Ricker_sim_tx <- Ricker_sim_fxn(data_tx_site, data_in$nwis_07332622)
Ricker_sim_pa <- Ricker_sim_fxn(data_pa_site, data_in$nwis_03025500)

# Save simulations.
saveRDS(Ricker_sim_mo, "data_working/Sim_MOriv_Ricker_jasm_2022_05_02.rds")
saveRDS(Ricker_sim_tx, "data_working/Sim_TXriv_Ricker_jasm_2022_05_02.rds")
saveRDS(Ricker_sim_pa, "data_working/Sim_PAriv_Ricker_jasm_2022_05_02.rds")

# And for each day at each site, I would like to calculate
# - mean GPP
# - 97.5% and 2.5% percentiles

# Going to pull out just the first site
data_mo_gpp <- Ricker_sim_mo[[1]]

# Calculate median and confidence intervals
median_gpp1 <- apply(data_mo_gpp, 1, median)
lowerci_gpp1 <- apply(data_mo_gpp, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp1 <- apply(data_mo_gpp, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp <- data_in$nwis_07061270$GPP

# Pull out original dates used
date <- data_in$nwis_07061270$date

# Bind into a single dataframe
df_mo_pred <- as.data.frame(cbind(median_gpp1, lowerci_gpp1, upperci_gpp1))

df_pred1 <- df_mo_pred %>%
  mutate(date = ymd(date),
         orig_gpp = orig_gpp)

# And plot
(gpp_plot1 <- ggplot(df_pred1 %>%
                       filter(date > "2014-12-31"), aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp1), color = "darkolivegreen2", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date") +
    geom_ribbon(aes(ymin = lowerci_gpp1,
                    ymax = upperci_gpp1),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12))) # 2015-2016?

# And now to pull out just the second site
data_tx_gpp <- Ricker_sim_tx[[1]]

# Calculate median and confidence intervals
median_gpp2 <- apply(data_tx_gpp, 1, median)
lowerci_gpp2 <- apply(data_tx_gpp, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp2 <- apply(data_tx_gpp, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp2 <- data_in$nwis_07332622$GPP

# Pull out original dates used
date2 <- data_in$nwis_07332622$date

# Bind into a single dataframe
df_tx_pred <- as.data.frame(cbind(median_gpp2, lowerci_gpp2, upperci_gpp2))

df_pred2 <- df_tx_pred %>%
  mutate(date = ymd(date2),
         orig_gpp = orig_gpp2)

# And plot
(gpp_plot2 <- ggplot(df_pred2 %>%
                       filter(date > "2015-12-31"), aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp2), color = "darkolivegreen2", size = 1.2) +
    labs(#y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date") +
    geom_ribbon(aes(ymin = lowerci_gpp2,
                    ymax = upperci_gpp2),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_blank())) # 2016-2017?

# And now to pull out just the third site
data_pa_gpp <- Ricker_sim_pa[[1]]

# Calculate median and confidence intervals
median_gpp3 <- apply(data_pa_gpp, 1, median)
lowerci_gpp3 <- apply(data_pa_gpp, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp3 <- apply(data_pa_gpp, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp3 <- data_in$nwis_03025500$GPP

# Pull out original dates used
date3 <- data_in$nwis_03025500$date

# Bind into a single dataframe
df_pa_pred <- as.data.frame(cbind(median_gpp3, lowerci_gpp3, upperci_gpp3))

df_pred3 <- df_pa_pred %>%
  mutate(date = ymd(date3),
         orig_gpp = orig_gpp3)

# And plot
(gpp_plot3 <- ggplot(df_pred3 %>%
                       filter(date > "2015-12-31"), aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp3), color = "darkolivegreen2", size = 1.2) +
    labs(#y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date") +
    geom_ribbon(aes(ymin = lowerci_gpp3,
                    ymax = upperci_gpp3),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_blank())) # 2016-2017

gpp_pred_fig <- gpp_plot1 | gpp_plot2 | gpp_plot3

# ggsave(("figures/teton_moresites/GPP_and_predGPP_3site.png"),
#        width = 40,
#        height = 10,
#        units = "cm"
# )

# End of script.
