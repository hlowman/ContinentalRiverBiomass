## Recovery of Stream Productivity following Disturbance
## Originally created: July 10, 2023
## Heili Lowman

#### READ ME ####

# The following set of scripts will walk through the steps necessary to
# prep and send data to Beartooth as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_XX" folders have been ignored using git.ignore.
# Links to the raw data sets are provided in the 01_NWIS_RiverSelection.R file.

# If you are accessing the code via GitHub, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running the 
# code below.

# This file will address remaining analyses to support responses to reviewer
# comments.

# The following script was run using the Pinyon server at the University of 
# Nevada Reno for speed.

#### Setup ####

# Load necessary packages.
lapply(c("tidybayes", "brms", "tidyverse", "lubridate", 
         "data.table", "GGally",
         "multcomp", "patchwork", "bayesplot",
         "modelsummary", "here", "nlme","loo"), 
       require, character.only=T)

#### Data ####

# Import necessary datasets.

# First, the model outputs.
dat_out <- readRDS("data_working/beartooth_181rivers_model_params_all_iterations_050823.rds")

# And, the original data used in modeling.
dat_in <- readRDS("data_working/df_181sites_Qmaxnorm_SavoySL.rds")

# And finally a dataset with long names for IDing purposes.
dat_names <- readRDS("data_working/NWIS_Info_181riv_HUC2_df_050923.rds") %>%
  dplyr::select(site_no, station_nm)

#### Reviewer 1 ####

##### Normalizing instead to Q10 #####

##### Pred GPP in the South Fork Iowa ######

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
  # and then predict forward for a site's timeseries (i.e., length(df$GPP))
  for(i in 1:length(pars$r)){
    simmat[,i] <- PM_Ricker(pars$r[i], pars$lambda[i], pars$s[i], pars$c[i], pars$sig_p[i], pars$sig_p[i], df)
    rmsemat[i] <- sqrt(sum((simmat[,i] - df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat, rmsemat)
  return(l)
  
}

# Adding the nRMSE calculation into the function above didn't play nicely with
# the list that existed, so calculating outside instead.
nRMSE_fxn <- function(df, df_orig){
  
  # Calculate the mean RMSE value for each site.
  nRMSE <- mean(df)/(max(df_orig$GPP) - min(df_orig$GPP))
  
}

# Take list containing all iterations of parameters and make into a df.
dat_out_df <- map_df(dat_out, ~as.data.frame(.x), .id="site_name")

# Trimming input and output datasets for the site of interest.
dat_out1df <- dat_out_df %>%
  filter(site_name %in% c("nwis_05451210"))

dat_in1df <- dat_in %>%
  filter(site_name %in% c("nwis_05451210"))

dat_in1 <- split(dat_in1df, dat_in1df$site_name)

# Re-simulating using all output iterations. Started ~1:28, Ended ~4:05
Ricker_sim1 <- Ricker_sim_fxn(dat_out1df, dat_in1df)

# And for each day, I would like to calculate 2.5%, 50%, and 97.5%tiles.

# Going to pull out just the predicted GPP values.
# So, making a list of odd numbers to pull out predGPP values.
data_1site_gpp <- Ricker_sim1[1]

# Calculate median and confidence intervals
quantile25 <- function(x){quantile(x, probs = 0.025, na.rm = TRUE)}
quantile975 <- function(x){quantile(x, probs = 0.975, na.rm = TRUE)}

pred_gpp1 <- lapply(data_1site_gpp, 
                    function(x) cbind(apply(x, 1, median),
                                      apply(x, 1, quantile25),
                                      apply(x, 1, quantile975)))

# Pull out original GPP values used and sequence #s (for plotting)
orig_gpp_date1 <- lapply(dat_in1, function(x) x %>% dplyr::select(date, GPP, seq,
                                                                  Q, Q_rel))

# Add names to confidence interval lists
my_names <- c("nwis_05451210")

names(pred_gpp1) <- my_names

pred_gpp1 <- lapply(pred_gpp1, function(x) as.data.frame(x) %>% 
                      rename("Median" = "V1",
                             "q2.5" = "V2",
                             "q97.5" = "V3")) # OMG YAY!!!!

# Bind into a single dataframe
keys <- unique(c(names(orig_gpp_date1), names(pred_gpp1)))
df_pred1 <- setNames(Map(cbind, orig_gpp_date1[keys], pred_gpp1[keys]), keys)

# And finally, calculate the normalized RMSE.
rmse1 <- Ricker_sim1[2]

nRMSE_1site <- mapply(nRMSE_fxn, rmse1, dat_in1)

# And plot
# South Fork Iowa River, IA
# First panel GPP
(gpp_plot_rmse1a <- ggplot(df_pred1$nwis_05451210 %>%
                             filter(date < ymd(as.character("2008-12-31"))), 
                           aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), color = "#609048", linewidth = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
        x = "Date") +
    scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2008-09-01"), y = 17.5,
             label = paste("nRMSE = ", round(nRMSE_1site[1], 
                           digits = 2)), size = 4) +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# First panel Q
(gpp_plot_rmse1a_q <- ggplot(df_pred1$nwis_05451210 %>%
                             filter(date < ymd(as.character("2008-12-31"))), 
                           aes(date, Q)) +
    geom_line(color = "#5792CC", linewidth = 1.2) +
    labs(y = expression('Discharge ('*m^3~s^-1*')'),
         x = "Date") +
    scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    ylim(0, 170) +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# Second panel
(gpp_plot_rmse1b <- ggplot(df_pred1$nwis_05451210 %>%
                             filter(date > ymd(as.character("2009-12-31"))), 
                           aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), color = "#609048", linewidth = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date") +
    scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# Second panel Q
(gpp_plot_rmse1b_q <- ggplot(df_pred1$nwis_05451210 %>%
                               filter(date > ymd(as.character("2009-12-31"))), 
                             aes(date, Q)) +
    geom_line(color = "#5792CC", linewidth = 1.2) +
    labs(y = expression('Discharge ('*m^3~s^-1*')'),
         x = "Date") +
    scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    ylim(0, 170) +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# Now, let's combine the above using patchwork.
(fig_SFIR <- (gpp_plot_rmse1a + gpp_plot_rmse1b +
                gpp_plot_rmse1a_q + gpp_plot_rmse1b_q) +
    plot_layout(nrow = 2) +
    plot_annotation(tag_levels = 'A'))

# And export.
# ggsave(fig_SFIR,
#        filename = "figures/beartooth_spring23/predGPP_SForkIowa_071023.jpg",
#        width = 30,
#        height = 15,
#        units = "cm")

##### Summary figures of sites lost during filtering #####

# Note, the revised version of the Supplementary figure displaying predicted 
# GPP at multiple sites can be found in the "14_Appendix_figures.R" script.

#### Reviewer 2 ####



# End of script.
