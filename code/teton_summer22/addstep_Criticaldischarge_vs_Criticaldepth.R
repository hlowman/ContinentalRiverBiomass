## Resilience of Stream Productivity to Disturbance
## August 31, 2022
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

# NOTE: This script is specifically designed to compared the Ricker model's
# calculations of "c" to critical discharge calculated based on
# geomorphic data provided by Jud Harvey.

#### Setup ####

## Load packages
lapply(c("tidyverse", "lubridate", "here", 
         "data.table", "pipeR", "patchwork",
         "ggplot2"), require, character.only=T)

#### Data Import ####

# Load site-level info (hypoxia and Appling datasets).
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
site <- fread("data_raw/site_data.tsv")

# Load site-level parameters dataset with all iterations included.
dat_out_q10 <- readRDS("data_working/teton_190rivers_model_site_params_all_iterations_082422.rds")
dat_out_max <- readRDS("data_working/teton_4rivers_model_site_params_all_iterations_090822.rds")
dat_out_max2 <- readRDS("data_working/teton_2rivers_model_site_params_all_iterations_091222.rds")

# Load in original data fed into the model.
dat_in <- readRDS("data_working/df_190sites_10yrQnorm_allSL.rds")
dat_in2list <- readRDS("data_working/list_4sites_meanmaxQnorm_allSL.rds")
dat_in2 <- map_df(dat_in2list, ~as.data.frame(.x), .id="site_name")
dat_in3list <- readRDS("data_working/list_2sites_meanmaxQnorm_allSL.rds")
dat_in3 <- map_df(dat_in3list, ~as.data.frame(.x), .id="site_name")

# Load in original model outputs to allow RMSE function to pull out each iteration set
# of parameters on its own.
# This data was part of a re-run where Q was normalized to Qmax instead of Q10yr flood.
data_out5 <- readRDS("data_teton/teton_4rivers_output_2022_09_08.rds")
data_out6 <- readRDS("data_teton/teton_2rivers_output_2022_09_12.rds")

#### Data Selection ####

# Geomorphic data is only available for the following sites:
# nwis_05451210: South Fork Iowa River, New Providence, IA
# nwis_05524500: Iroquois River, Foresman, IN
# nwis_05579630: Kickapoo Creek, Bloomington, IL
# nwis_01645704: Difficult Run, Fairfax, VA
# nwis_0165389205: Accotink Creek, Ranger Road, VA
# nwis_05515500: Kankakee River, Davis, IN

# Filtering by 4 sites
my_list <- c("nwis_05451210","nwis_05524500",
             "nwis_05579630","nwis_01645704")

# Filtering by all 6 sites
my_list6 <- c("nwis_05451210","nwis_05524500",
              "nwis_05579630","nwis_01645704",
              "nwis_0165389205", "nwis_05515500")

dat_out_q10_6 <- dat_out_q10 %>%
  filter(site_name %in% my_list6)

# And calculate normalized c statistic values
c_summary_q10 <- dat_out_q10_6 %>%
  group_by(site_name) %>%
  summarize(mean_c = mean(c),
            med_c = median(c),
            q2.5_c = quantile(c, c(.025)),
            q97.5_c = quantile(c, c(.975))) %>%
  ungroup()

# Identify ten year flood values to back calculate c in cfs.
Q10_6 <- dat_in %>%
  filter(site_name %in% my_list6) %>%
  distinct(site_name,RI_10yr_Q,RI_10yr_Q_cms)

# Join these datasets.
c_Q10_6 <- left_join(c_summary_q10, Q10_6) # joined by "site_name"

# Calculate "c" in cubic meters per second.
c_Q10_6 <- c_Q10_6 %>%
  mutate(mean_c_cms = mean_c * RI_10yr_Q_cms,
         q50_c_cms = med_c * RI_10yr_Q_cms,
         q2.5_c_cms = q2.5_c * RI_10yr_Q_cms,
         q97.5_c_cms = q97.5_c * RI_10yr_Q_cms)

# Convert to cfs.
c_Q10_6 <- c_Q10_6 %>%
  mutate(mean_c_cfs = mean_c_cms * 35.3147,
         q50_c_cfs = q50_c_cms * 35.3147,
         q2.5_c_cfs = q2.5_c_cms * 35.3147,
         q97.5_c_cfs = q97.5_c_cms * 35.3147)

# Trim and export data.
c_Q10_6_trim <- c_Q10_6 %>%
  dplyr::select(site_name, mean_c_cfs, q2.5_c_cfs, q50_c_cfs, q97.5_c_cfs)

# write_csv(c_Q10_6_trim, 
#           "data_working/critical_discharge_cfs_6_sites_091222.csv")

#### Data Prep ####

# Re-running these 6 sites and normalizing by MAX Q rather than 10yr flood Q.
# Tried average Q, but the model freaked out because values were not between 0 and 1.

# First, filter available data by the sites of interest.
dat_in_4 <- dat_in %>%
  filter(site_name %in% my_list)

# And calculate mean/max Q for time periods for which we have data.
q_stats <- dat_in_4 %>%
  group_by(site_name) %>%
  summarize(meanQ = mean(Q),
            maxQ = max(Q))

# Join this info to the main dataset.
dat_in_4 <- left_join(dat_in_4, q_stats)

# And create the new columns with Q normalized to the mean & max.
dat_in_4 <- dat_in_4 %>%
  mutate(Q_rel_mean = Q/meanQ,
         Q_rel_max = Q/maxQ)

# Make as list and export for use on Teton.
dat4_list <- split(dat_in_4, dat_in_4$.id)

#saveRDS(dat4_list, "data_working/list_4sites_meanmaxQnorm_allSL.rds")

# And do the same for 2 additional sites per Jud's request.
# 0165389205 - Accotink Creek near Ranger Road, VA
# 05515500 - Kankakee River at Davis, IN

my_list2 <- c("nwis_0165389205", "nwis_05515500")

# First, filter available data by the 4 sites of interest.
dat_in_2 <- dat_in %>%
  filter(site_name %in% my_list2)

# And calculate mean Q for time periods for which we have data.
q_stats2 <- dat_in_2 %>%
  group_by(site_name) %>%
  summarize(meanQ = mean(Q),
            maxQ = max(Q))

# Join this info to the main dataset.
dat_in_2 <- left_join(dat_in_2, q_stats2)

# And create the new column with Q normalized to the mean.
dat_in_2 <- dat_in_2 %>%
  mutate(Q_rel_mean = Q/meanQ,
         Q_rel_max = Q/maxQ)

# Make as list and export for use on Teton.
dat2_list <- split(dat_in_2, dat_in_2$.id)

#saveRDS(dat2_list, "data_working/list_2sites_meanmaxQnorm_allSL.rds")

#### Data Processing ####

# Calculate normalized c statistic values.
# First, bind together data_out_max (4 sites) and data_out_max2 (2 sites).
dat_out_max_6 <- rbind(dat_out_max, dat_out_max2)

c_summary_max <- dat_out_max_6 %>%
  group_by(site_name) %>%
  summarize(mean_c = mean(c),
            med_c = median(c),
            q2.5_c = quantile(c, c(.025)),
            q97.5_c = quantile(c, c(.975))) %>%
  ungroup()

# Identify maximum discharge values to back calculate c in cfs.
Qmax_4 <- dat_in2 %>%
  filter(site_name %in% my_list) %>%
  distinct(site_name,maxQ)

Qmax_2 <- dat_in3 %>%
  filter(site_name %in% my_list2) %>%
  distinct(site_name,maxQ)

Qmax_6 <- rbind(Qmax_4, Qmax_2)

# Join these datasets.
c_Qmax_6 <- left_join(c_summary_max, Qmax_6) # joined by "site_name"

# Calculate "c" in cubic meters per second.
c_Qmax_6 <- c_Qmax_6 %>%
  mutate(mean_c_cms = mean_c * maxQ,
         q50_c_cms = med_c * maxQ,
         q2.5_c_cms = q2.5_c * maxQ,
         q97.5_c_cms = q97.5_c * maxQ)

# Convert to cfs.
c_Qmax_6 <- c_Qmax_6 %>%
  mutate(mean_c_cfs = mean_c_cms * 35.3147,
         q50_c_cfs = q50_c_cms * 35.3147,
         q2.5_c_cfs = q2.5_c_cms * 35.3147,
         q97.5_c_cfs = q97.5_c_cms * 35.3147)

# Trim and export data.
c_Qmax_6_trim <- c_Qmax_6 %>%
  dplyr::select(site_name, mean_c_cfs, q2.5_c_cfs, q50_c_cfs, q97.5_c_cfs)

# write_csv(c_Qmax_6_trim, 
#           "data_working/critical_discharge_cfs_qmax_6_sites_091222.csv")

# Plot for lab meeting next week.
c_Qmax_6_trim$Normalization <- "Qmax"
c_Q10_6_trim$Normalization <- "Q10"

c_all <- rbind(c_Qmax_6_trim, c_Q10_6_trim)

(fig_c <- ggplot(c_all, aes(x = mean_c_cfs, y = site_name, color = Normalization)) +
    geom_point(size = 5, alpha = 0.75, position = position_dodge(width = -0.5)) +
    geom_errorbarh(aes(xmin = q2.5_c_cfs, xmax = q97.5_c_cfs), 
                   height = 0.25,
                   position = position_dodge(width = -0.5)) +
    labs(y = "Site", x = "c (cfs)") +
    scale_color_manual(values = c("#7AC9B7","#3793EC")) +
    theme_bw())

# export exploratory figures
# ggsave(("figures/teton_summer22/c_q10_vs_qmaxnorm_fig.png"),
#        width = 16,
#        height = 10,
#        units = "cm"
# )

#### RMSE/GPP Prediction ####

# Calculate RMSE of predicted GPP using qmax normalized 6 site parameter estimates.

# Previously, I was just pulling model diagnostics. The below workflow should now
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
  tQ <- df$Q_rel_max # discharge standardized to max value
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
  
  # extract parameters from STAN output
  pars <- extract(output, c("r", "lambda", "s", "c", "B", "P", "pred_GPP", "sig_p", "sig_o"))
  #pars <- output
  
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

# And finally, apply the function to my 6 sites.
Ricker_sim1 <- Ricker_sim_fxn(data_out5$nwis_01645704, dat_in2list$nwis_01645704)
Ricker_sim2 <- Ricker_sim_fxn(data_out5$nwis_05451210, dat_in2list$nwis_05451210)
Ricker_sim3 <- Ricker_sim_fxn(data_out5$nwis_05524500, dat_in2list$nwis_05524500)
Ricker_sim4 <- Ricker_sim_fxn(data_out5$nwis_05579630, dat_in2list$nwis_05579630)
Ricker_sim5 <- Ricker_sim_fxn(data_out6$nwis_0165389205, dat_in3list$nwis_0165389205)
Ricker_sim6 <- Ricker_sim_fxn(data_out6$nwis_05515500, dat_in3list$nwis_05515500)

# Save simulations.
# saveRDS(Ricker_sim1, "data_working/Sim_Ricker_01645704_2022_09_26.rds")
# saveRDS(Ricker_sim2, "data_working/Sim_Ricker_05451210_2022_09_26.rds")
# saveRDS(Ricker_sim3, "data_working/Sim_Ricker_05524500_2022_09_26.rds")
# saveRDS(Ricker_sim4, "data_working/Sim_Ricker_05579630_2022_09_26.rds")
# saveRDS(Ricker_sim5, "data_working/Sim_Ricker_0165389205_2022_09_26.rds")
# saveRDS(Ricker_sim6, "data_working/Sim_Ricker_05515500_2022_09_26.rds")

# And for each day at each site, I would like to calculate
# - mean GPP
# - 97.5% and 2.5% percentiles

#### Difficult Run, VA ####

# Going to pull out just the first site
data_sim1_gpp <- Ricker_sim1[[1]]
data_sim1_rmse <- Ricker_sim1[[2]]

# Calculate median and confidence intervals
median_gpp1 <- apply(data_sim1_gpp, 1, median)
lowerci_gpp1 <- apply(data_sim1_gpp, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp1 <- apply(data_sim1_gpp, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp1 <- dat_in2list$nwis_01645704$GPP

# Pull out original dates used
date1 <- dat_in2list$nwis_01645704$date

# Pull out original discharge values used
discharge1 <- dat_in2list$nwis_01645704$Q

# Bind into a single dataframe
df_sim1_pred <- as.data.frame(cbind(median_gpp1, lowerci_gpp1, upperci_gpp1))

df_pred1 <- df_sim1_pred %>%
  mutate(date = ymd(date1),
         orig_gpp = orig_gpp1,
         discharge = discharge1)

# calculate RMSE
mean_sim1_rmse <- mean(data_sim1_rmse)

# calculate normalized RMSE
nrmse1 <- mean_sim1_rmse/(max(df_pred1$orig_gpp)-min(df_pred1$orig_gpp))

# And plot
(gpp_plot1 <- ggplot(df_pred1, aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp1), color = "darkolivegreen2", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Difficult Run, Fairfax, VA\nnwis 01645704") +
    geom_ribbon(aes(ymin = lowerci_gpp1,
                    ymax = upperci_gpp1),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    annotate(geom = "text", x = date("2016-07-01"), y = 10, 
             label = paste("Normalized RMSE = ", 
                           round(nrmse1, digits = 2), 
                           sep = "")) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)))

# Plot hydrograph as well.
(gpp_plot1.2 <- ggplot(df_pred1, aes(date, discharge1)) +
    geom_line(aes(date, discharge1), color = "deepskyblue4", size = 1.2) +
    labs(y = "Discharge (cms)",
         x = "Date") +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)))

(VA_plot <- gpp_plot1 / gpp_plot1.2)

# Export to send to Jud.
# ggsave(("figures/teton_summer22/GPP_and_predGPP_difficultrun.png"),
#        width = 30,
#        height = 20,
#        units = "cm"
# )

#### South Fork, IA ####

# Going to pull out just the second site
data_sim2_gpp <- Ricker_sim2[[1]]
data_sim2_rmse <- Ricker_sim2[[2]]

# Calculate median and confidence intervals
median_gpp2 <- apply(data_sim2_gpp, 1, median)
lowerci_gpp2 <- apply(data_sim2_gpp, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp2 <- apply(data_sim2_gpp, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp2 <- dat_in2list$nwis_05451210$GPP

# Pull out original dates used
date2 <- dat_in2list$nwis_05451210$date

# Pull out original discharge values used
discharge2 <- dat_in2list$nwis_05451210$Q

# Bind into a single dataframe
df_sim2_pred <- as.data.frame(cbind(median_gpp2, lowerci_gpp2, upperci_gpp2))

df_pred2 <- df_sim2_pred %>%
  mutate(date = ymd(date2),
         orig_gpp = orig_gpp2,
         discharge = discharge2)

# calculate RMSE
mean_sim2_rmse <- mean(data_sim2_rmse)

# calculate normalized RMSE
nrmse2 <- mean_sim2_rmse/(max(df_pred2$orig_gpp)-min(df_pred2$orig_gpp))

# And plot
(gpp_plot2 <- ggplot(df_pred2, aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp2), color = "darkolivegreen2", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "South Fork Iowa River, New Providence, IA\nnwis 05451210") +
    geom_ribbon(aes(ymin = lowerci_gpp2,
                    ymax = upperci_gpp2),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    annotate(geom = "text", x = date("2010-09-01"), y = 20, 
             label = paste("Normalized RMSE = ", 
                           round(nrmse2, digits = 2), 
                           sep = "")) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)))

# Plot hydrograph as well.
(gpp_plot2.2 <- ggplot(df_pred2, aes(date, discharge)) +
    geom_line(aes(date, discharge), color = "deepskyblue4", size = 1.2) +
    labs(y = "Discharge (cms)",
         x = "Date") +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)))

(IA_plot <- gpp_plot2 / gpp_plot2.2)

# Export to send to Jud.
# ggsave(("figures/teton_summer22/GPP_and_predGPP_southfork.png"),
#        width = 30,
#        height = 20,
#        units = "cm"
# )

#### Iroquois River, IN ####

# Going to pull out just the third site
data_sim3_gpp <- Ricker_sim3[[1]]
data_sim3_rmse <- Ricker_sim3[[2]]

# Calculate median and confidence intervals
median_gpp3 <- apply(data_sim3_gpp, 1, median)
lowerci_gpp3 <- apply(data_sim3_gpp, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp3 <- apply(data_sim3_gpp, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp3 <- dat_in2list$nwis_05524500$GPP

# Pull out original dates used
date3 <- dat_in2list$nwis_05524500$date

# Bind into a single dataframe
df_sim3_pred <- as.data.frame(cbind(median_gpp3, lowerci_gpp3, upperci_gpp3))

df_pred3 <- df_sim3_pred %>%
  mutate(date = ymd(date3),
         orig_gpp = orig_gpp3)

# calculate RMSE
mean_sim3_rmse <- mean(data_sim3_rmse)

# And plot
(gpp_plot3 <- ggplot(df_pred3, aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp3), color = "darkolivegreen2", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Iroquois River, Foresman, IN") +
    geom_ribbon(aes(ymin = lowerci_gpp3,
                    ymax = upperci_gpp3),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    annotate(geom = "text", x = date("2016-11-01"), y = 10, 
             label = paste("RMSE = ", round(mean_sim3_rmse, digits = 2), sep = "")) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)))

#### Kickapoo Creek, IL ####

# Going to pull out just the fourth site
data_sim4_gpp <- Ricker_sim4[[1]]
data_sim4_rmse <- Ricker_sim4[[2]]

# Calculate median and confidence intervals
median_gpp4 <- apply(data_sim4_gpp, 1, median)
lowerci_gpp4 <- apply(data_sim4_gpp, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp4 <- apply(data_sim4_gpp, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp4 <- dat_in2list$nwis_05579630$GPP

# Pull out original dates used
date4 <- dat_in2list$nwis_05579630$date

# Bind into a single dataframe
df_sim4_pred <- as.data.frame(cbind(median_gpp4, lowerci_gpp4, upperci_gpp4))

df_pred4 <- df_sim4_pred %>%
  mutate(date = ymd(date4),
         orig_gpp = orig_gpp4)

# calculate RMSE
mean_sim4_rmse <- mean(data_sim4_rmse)

# And plot
(gpp_plot4 <- ggplot(df_pred4, aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp4), color = "darkolivegreen2", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Kickapoo Creek, Bloomington, IL") +
    geom_ribbon(aes(ymin = lowerci_gpp4,
                    ymax = upperci_gpp4),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    annotate(geom = "text", x = date("2012-07-01"), y = 14, 
             label = paste("RMSE = ", round(mean_sim4_rmse, digits = 2), sep = "")) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)))

#### Accotink Creek, VA ####
# Going to pull out just the fifth site
data_sim5_gpp <- Ricker_sim5[[1]]
data_sim5_rmse <- Ricker_sim5[[2]]

# Calculate median and confidence intervals
median_gpp5 <- apply(data_sim5_gpp, 1, median)
lowerci_gpp5 <- apply(data_sim5_gpp, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp5 <- apply(data_sim5_gpp, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp5 <- dat_in3list$nwis_0165389205$GPP

# Pull out original dates used
date5 <- dat_in3list$nwis_0165389205$date

# Pull out original discharge used
discharge5 <- dat_in3list$nwis_0165389205$Q

# Bind into a single dataframe
df_sim5_pred <- as.data.frame(cbind(median_gpp5, lowerci_gpp5, upperci_gpp5))

df_pred5 <- df_sim5_pred %>%
  mutate(date = ymd(date5),
         orig_gpp = orig_gpp5,
         discharge = discharge5)

# calculate RMSE
mean_sim5_rmse <- mean(data_sim5_rmse)

# calculate normalized RMSE
nrmse5 <- mean_sim5_rmse/(max(df_pred5$orig_gpp)-min(df_pred5$orig_gpp))

# And plot
(gpp_plot5 <- ggplot(df_pred5, aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp5), color = "darkolivegreen2", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Accotink Creek, Ranger Road, VA\nnwis 0165389205") +
    geom_ribbon(aes(ymin = lowerci_gpp5,
                    ymax = upperci_gpp5),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    annotate(geom = "text", x = date("2013-12-01"), y = 6.5, 
             label = paste("Normalized RMSE = ", 
                           round(nrmse5, digits = 2), 
                           sep = "")) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)))

# Plot hydrograph as well.
(gpp_plot5.2 <- ggplot(df_pred5, aes(date, discharge)) +
    geom_line(aes(date, discharge), color = "deepskyblue4", size = 1.2) +
    labs(y = "Discharge (cms)",
         x = "Date") +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)))

(VA_plot2 <- gpp_plot5 / gpp_plot5.2)

# Export to send to Jud.
# ggsave(("figures/teton_summer22/GPP_and_predGPP_accotink.png"),
#        width = 30,
#        height = 20,
#        units = "cm"
# )

#### Kankakee River, IN ####

# Going to pull out just the sixth site
data_sim6_gpp <- Ricker_sim6[[1]]
data_sim6_rmse <- Ricker_sim6[[2]]

# Calculate median and confidence intervals
median_gpp6 <- apply(data_sim6_gpp, 1, median)
lowerci_gpp6 <- apply(data_sim6_gpp, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp6 <- apply(data_sim6_gpp, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp6 <- dat_in3list$nwis_05515500$GPP

# Pull out original dates used
date6 <- dat_in3list$nwis_05515500$date

# Bind into a single dataframe
df_sim6_pred <- as.data.frame(cbind(median_gpp6, lowerci_gpp6, upperci_gpp6))

df_pred6 <- df_sim6_pred %>%
  mutate(date = ymd(date6),
         orig_gpp = orig_gpp6)

# calculate RMSE
mean_sim6_rmse <- mean(data_sim6_rmse)

# And plot
(gpp_plot6 <- ggplot(df_pred6, aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp6), color = "darkolivegreen2", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Kankakee River, Davis, IN") +
    geom_ribbon(aes(ymin = lowerci_gpp6,
                    ymax = upperci_gpp6),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    annotate(geom = "text", x = date("2016-12-01"), y = 3.25, 
             label = paste("RMSE = ", round(mean_sim6_rmse, digits = 2), sep = "")) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)))

(gpp_pred_fig <- gpp_plot1 / gpp_plot2 / gpp_plot3 / gpp_plot4 / gpp_plot5 / gpp_plot6)

# ggsave(("figures/teton_summer22/GPP_and_predGPP_6site.png"),
#        width = 40,
#        height = 60,
#        units = "cm"
# )

# End of script.
