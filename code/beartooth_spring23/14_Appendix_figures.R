## Recovery of Stream Productivity following Disturbance
## Originally created: May 12, 2023
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

#### Setup ####

# Load necessary packages.
lapply(c("tidybayes", "brms", "tidyverse", "lubridate", 
         "data.table", "GGally", "patchwork", "here",
         "reshape2", "ggExtra", "ggbreak", "sf",
         "maps", "mapproj"), 
       require, character.only=T)

#### Data ####

# Import necessary datasets.

# First, the data for the maximum accrual (amax) models.
dat_amax <- readRDS("data_working/amax_covariates_152sites_070523.rds")

# Next, the data for the Qc:Q2yr models.
dat_Qc <- readRDS("data_working/Qc_covariates_138sites_070523.rds")

# And, the original data used in modeling.
dat_in <- readRDS("data_working/df_181sites_Qmaxnorm_SavoySL.rds")

# As well as the model outputs.
dat_out <- readRDS("data_working/beartooth_181rivers_model_params_all_iterations_050823.rds")

# And all covariate data at all sites.
dat_cov <- readRDS("data_working/covariate_data_181sites_070523.rds")

#### Timeseries Length Appendix Figure ####

# Making a histogram of the lengths of the timeseries for an appendix figure.
site_lengths <- dat_in %>%
  group_by(site_name) %>%
  summarize(days = n()) %>% # count the number of rows/days
  ungroup()

# Adding bins for more intuitive plotting.
site_lengths <- site_lengths %>%
  mutate(bin = factor(case_when(days < 180 ~ "3-6 months",
                                days >= 180 & days < 365 ~ "6 mos. - 1 year",
                                days >= 365 & days < 730 ~ "1-2 years",
                                days >= 730 & days < 1460 ~ "2-4 years",
                                days >= 1460 & days < 2190 ~ "4-6 years",
                                days >= 2190 & days < 2920 ~ "6-8 years",
                                days >= 2920 ~ "8-10 years"),
                      levels = c("3-6 months", "6 mos. - 1 year",
                                 "1-2 years", "2-4 years",
                                 "4-6 years", "6-8 years",
                                 "8-10 years")))

(ts_hist <- ggplot(site_lengths, aes(x = bin)) +
    geom_histogram(stat = "count", fill = "#7AC9B7") +
    labs(x = "Length of Timeseries",
         y = "Number of Sites") +
    theme_bw() +
    theme(text = element_text(size = 10)))

# And export.
# ggsave(ts_hist,
#        filename = "figures/beartooth_spring23/TS_length_hist_051223.jpg",
#        width = 14,
#        height = 7,
#        units = "cm")

# And calculate median timeseries length for inclusion in methods.
median(site_lengths$days) # 870

# And calculate site-years for inclusion in methods.
dat_sy <- dat_in %>%
  group_by(site_name, year) %>%
  summarize(days = n()) %>%
  ungroup() # 684

mean(dat_sy$days) # mean of 283 days/year

#### Site Map Appendix Figure ####

# make base US map
states <- map_data("state")
states_sf <- st_as_sf(states,
                      coords = c("long", "lat"),
                      remove = F,
                      crs = 4326)

# make data sf object
sites_sf <- st_as_sf(dat_cov,
                     coords = c("Lon_WGS84", 
                                "Lat_WGS84"), # always put lon (x) first
                     remove = F, # leave the lat/lon columns in too
                     crs = 4326) # projection: WGS84
# see for more mapping info: https://www.nceas.ucsb.edu/sites/default/files/2020-04/OverviewCoordinateReferenceSystems.pdf

# site map colored by model used
(sitemap <- ggplot(states_sf) + # base plot
    geom_polygon(aes(x = long, y = lat, group = group), 
                 fill = "white", color = "black") + # map of states
    geom_point(data = sites_sf, aes(x = Lon_WGS84, y = Lat_WGS84), 
               fill = "#7AC9B7", shape = 21, 
               size = 4, alpha = 0.75) + # map of sites
    theme_classic() + # remove grid
    labs(x = "Longitude",
         y = "Latitude") +
    coord_map(projection = "albers", lat0 = 39, lat1 = 45))

# Checked metadata for max/min lat/lon to be sure no states outside CONUS
# included in modeling.

ggsave(plot = sitemap,
       filename = "figures/beartooth_spring23/map_181sites.jpg",
       width = 20,
       height = 12,
       units = "cm")

#### Accrual Appendix Figure ####

# MAY vs. GPP:
(figa1 <- ggplot(dat_amax, aes(x = meanGPP, y = yield_med)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#96AA8B") +
    geom_linerange(alpha = 0.8, 
                   color = "#96AA8B",
                   aes(ymin = yield_2.5_ed, ymax = yield_97.5)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = expression(a[max]),
         x = expression(Mean~Daily~GPP~(gO[2]~m^-2~d^-1))) +
    theme_bw())

# MAY vs. cvQ:
(figa2 <- ggplot(dat_amax, aes(x = cvQ, y = yield_med)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#8A9C7E") +
    geom_linerange(alpha = 0.8, 
                   color = "#8A9C7E",
                   aes(ymin = yield_2.5_ed, ymax = yield_97.5)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(CV[Q]),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. daily light:
(figa3 <- ggplot(dat_amax, aes(x = meanLight, y = yield_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#7B8B6E") +
    geom_linerange(alpha = 0.8, 
                   color = "#7B8B6E",
                   aes(ymin = yield_2.5_ed, ymax = yield_97.5)) +
    scale_y_log10() +
    labs(x = expression(Mean~Daily~PAR~(mol~m^-2~d^-1)),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. annual exceedances:
(figa4 <- ggplot(dat_amax, aes(x = exc_y, y = yield_med)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#6C7A5D") +
    geom_linerange(alpha = 0.8, 
                   color = "#6C7A5D",
                   aes(ymin = yield_2.5_ed, ymax = yield_97.5)) +
    scale_y_log10() +
    labs(y = expression(a[max]),
         x = expression(Mean~Annual~Exceedances~of~Q[c])) +
    theme_bw())

# MAY vs. HUC: 
(figa5 <- ggplot(dat_amax, aes(x = huc_2, y = yield_med)) +
    geom_boxplot(alpha = 0.6, color = "black", fill = "#6C7A5D") +
    scale_y_log10() +
    labs(x = expression(Regional~Hydrological~Unit~Code),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. summer temp:
(figa6 <- ggplot(dat_amax, aes(x = summermeanTemp, y = yield_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#5C694C") +
    geom_linerange(alpha = 0.8, 
                   color = "#5C694C",
                   aes(ymin = yield_2.5_ed, ymax = yield_97.5)) +
    scale_y_log10() +
    labs(x = paste0("Mean Summer Temperature (", '\u00B0', "C)"),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. dams:
(figa7 <- ggplot(dat_amax %>%
                    # Creating the new categorical dam column we modeled by.
                    mutate(Dam_binary = factor(case_when(
                      Dam %in% c("50", "80", "95") ~ "0", # Potential
                      Dam == "0" ~ "1", # Certain
                      TRUE ~ NA))) %>%
                    drop_na(Dam_binary), 
                   aes(x = Dam_binary, y = yield_med)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#4D583C", color = "black") +
    scale_y_log10() +
    scale_x_discrete(labels = c("5-50%", "100%")) +
    labs(x = expression(Likelihood~of~Dam~Influence),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. stream width:
(figa8 <- ggplot(dat_amax, aes(x = width_med, y = yield_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#454F34") +
    geom_linerange(alpha = 0.8, 
                   color = "#454F34",
                   aes(ymin = yield_2.5_ed, ymax = yield_97.5)) +
    scale_x_log10() + 
    scale_y_log10() +
    labs(x = expression(Stream~Width~(m)),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. Road density:
(figa9 <- ggplot(dat_amax, 
                  aes(x = NHD_RdDensWs, y = yield_med)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#3A422D") +
    geom_linerange(alpha = 0.8, 
                   color = "#3A422D",
                   aes(ymin = yield_2.5_ed, ymax = yield_97.5)) +
    scale_y_log10() +
    labs(x = expression(Road~Density~(km/km^2)),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. NO3:
(figa10 <- ggplot(dat_amax, aes(x = Nitrate, y = yield_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#343B29") +
    geom_linerange(alpha = 0.8, 
                   color = "#343B29",
                   aes(ymin = yield_2.5_ed, ymax = yield_97.5)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Mean~Nitrate~(mg/L~NO[3]-N)),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. P:
(figa11 <- ggplot(dat_amax, aes(x = Phosphorus, y = yield_med)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#2F3525") +
    geom_linerange(alpha = 0.8, 
                   color = "#2F3525",
                   aes(ymin = yield_2.5_ed, ymax = yield_97.5)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Mean~Dissolved~Phosphorus~(mg/L~P)),
         y = expression(a[max])) +
    theme_bw())

# Combine figures above.
(fig_yield_med <- figa1 + figa2 + figa3 +
    figa4 + figa5 + figa6 + 
    figa7 + figa8 + figa9 +
    figa10 + figa11 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 4))

# And export for use in the appendix.
# ggsave(fig_yield_med,
#        filename = "figures/beartooth_spring23/amax_11panel_051223.jpg",
#        width = 22,
#        height = 26,
#        units = "cm") # n = 152

#### Qc Appendix Figure ####

# Create columns necessary for plotting.
dat_Qc <- dat_Qc %>%
  mutate(Qc_Q2yr = Qc/RI_2yr_Q_cms,
         Qc_Q2yr2.5 = Qc2.5/RI_2yr_Q_cms,
         Qc_Q2yr97.5 = Qc97.5/RI_2yr_Q_cms)

# Mean daily GPP vs. Qc:Q2: X axis LOG SCALED
(figq1 <- ggplot(dat_Qc, aes(x = meanGPP, y = Qc_Q2yr)) +
   geom_point(alpha = 0.8, size = 3,
              color = "#486999") +
   geom_linerange(alpha = 0.8, 
                  color = "#486999",
                  aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
   geom_hline(yintercept = 1, linetype = "dashed") +
   scale_x_log10() +
   scale_y_log10() + 
   labs(y = expression(Q[c]:Q[2~yr]),
        x = expression(Mean~Daily~GPP~(gO[2]~m^-2~d^-1))) +
   theme_bw())

# CV of Discharge vs. Qc:Q2: note, x axis on LOG SCALE
(figq2 <- ggplot(dat_Qc, aes(x = cvQ, y = Qc_Q2yr)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#405F8A") +
    geom_linerange(alpha = 0.8, 
                   color = "#405F8A",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    scale_y_log10() + 
    labs(x = expression(CV[Q]),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Mean annual exceedances vs. Qc:Q2
(figq3 <- ggplot(dat_Qc, aes(x = exc_y, y = Qc_Q2yr)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#38557A") +
    geom_linerange(alpha = 0.8, 
                   color = "#38557A",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    scale_y_log10() + 
    labs(x = expression(Mean~Annual~Exceedances~of~Q[c]),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# HUC vs. Qc:Q2: 
(figq4 <- ggplot(dat_Qc, aes(x = huc_2, y = Qc_Q2yr)) +
    geom_boxplot(alpha = 0.6, color = "black", fill = "#38557A") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() +
    labs(x = expression(Regional~Hydrological~Unit~Code),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Effect of Dams
(figq5 <- ggplot(dat_Qc %>%
                       # Creating the new categorical dam column we modeled by.
                       mutate(Dam_binary = factor(case_when(
                         Dam %in% c("50", "80", "95") ~ "0", # Potential
                         Dam == "0" ~ "1", # Certain
                         TRUE ~ NA))) %>%
                       drop_na(Dam_binary), 
                     aes(x = Dam_binary, y = Qc_Q2yr)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#273C57", color = "black") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() + 
    scale_x_discrete(labels = c("5-50%", "100%")) +
    labs(x = expression(Likelihood~of~Influence~by~Dams),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Stream Width vs. c: note, x axis LOG SCALED
(figq6 <- ggplot(dat_Qc, aes(x = width_med, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#1E2F46") +
    geom_linerange(alpha = 0.8, 
                   color = "#1E2F46",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_log10() + 
    scale_y_log10() + 
    labs(x = expression(Stream~Width~(m)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Road density in the watershed
(figq7 <- ggplot(dat_Qc, aes(x = NHD_RdDensWs, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#304969") +
    geom_linerange(alpha = 0.8, 
                   color = "#304969",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() + 
    labs(x = expression(Road~Density~(km/km^2)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Combine figures above and export for supplemental figure.
(fig_qcq2_supp <- figq1 + figq2 + figq3 +
    figq4 + figq5 + figq6 + figq7 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 3))

# ggsave(fig_qcq2_supp,
#        filename = "figures/beartooth_spring23/QcQ2_7panel_051223.jpg",
#        width = 22,
#        height = 21,
#        units = "cm") # n = 138

#### GPP Appendix Figure ####

# 8-paneled plot demonstrating model fit across a variety of site-types.

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

my8sites <- c("nwis_0166818623", "nwis_02217643",
              "nwis_06893350", "nwis_07075250",
              "nwis_05082500", "nwis_13013650",
              "nwis_04176500", "nwis_08374550")

# Take list containing all iterations of parameters and make into a df.
dat_out_df <- map_df(dat_out, ~as.data.frame(.x), .id="site_name")

# Trimming input and output datasets for the sites of interest.
dat_out8df <- dat_out_df %>%
  filter(site_name %in% my8sites)

dat_out8 <- split(dat_out8df, dat_out8df$site_name)

dat_in8df <- dat_in %>%
  filter(site_name %in% my8sites)

dat_in8 <- split(dat_in8df, dat_in8df$site_name)

# Re-simulating using all output iterations. Started ~3:45, Ended ~4:05
Ricker_sim8sites <- mapply(Ricker_sim_fxn, dat_out8, dat_in8)

# And for each day, I would like to calculate
# - median GPP
# - 97.5% and 2.5% percentiles

# Going to pull out just the predicted GPP values.
# So, making a list of odd numbers to pull out predGPP values (see above for reasoning).
my_values <- seq(from = 1, to = 16, by = 2)
data_8site_gpp <- Ricker_sim8sites[my_values]

# Calculate median and confidence intervals
quantile25 <- function(x){quantile(x, probs = 0.025, na.rm = TRUE)}
quantile975 <- function(x){quantile(x, probs = 0.975, na.rm = TRUE)}

pred_gpp8 <- lapply(data_8site_gpp, 
                    function(x) cbind(apply(x, 1, median),
                                      apply(x, 1, quantile25),
                                      apply(x, 1, quantile975)))

# Pull out original GPP values used and sequence #s (for plotting)
orig_gpp_date8 <- lapply(dat_in8, function(x) x %>% select(date, GPP, seq))

# Add names to confidence interval lists
my_names <- c("nwis_0166818623", "nwis_02217643", "nwis_04176500", "nwis_05082500", "nwis_06893350", "nwis_07075250", "nwis_08374550", "nwis_13013650")

names(pred_gpp8) <- my_names
pred_gpp8 <- lapply(pred_gpp8, function(x) as.data.frame(x) %>% 
                      rename("Median" = "V1",
                             "q2.5" = "V2",
                             "q97.5" = "V3")) # OMG YAY!!!!

# Bind into a single dataframe
keys <- unique(c(names(orig_gpp_date8), names(pred_gpp8)))
df_pred8 <- setNames(Map(cbind, orig_gpp_date8[keys], pred_gpp8[keys]), keys)

# And finally, calculate the normalized RMSE.
my_values2 <- seq(from = 2, to = 16, by = 2)
rmse8 <- Ricker_sim8sites[my_values2]

nRMSE_8site <- mapply(nRMSE_fxn, rmse8, dat_in8)

# And plot
# Mill Creek, VA
(gpp_plot8.1 <- ggplot(df_pred8$nwis_0166818623, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Small Order - Steady Flow - High Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2016-06-01"), y = 13,
             label = paste("nRMSE = ",round(nRMSE_8site[1], 
                                            digits = 2)), size = 4) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10), 
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# Parks Creek, GA
(gpp_plot8.2 <- ggplot(df_pred8$nwis_02217643, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Small Order - Steady Flow - Low Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2016-11-01"), y = 6,
             label = paste("nRMSE = ",round(nRMSE_8site[2], 
                                            digits = 2)), size = 4) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10), 
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# Tomahawk Creek, KS
(gpp_plot8.3 <- ggplot(df_pred8$nwis_06893350, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Small Order - Turbulent Flow - High Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2011-12-01"), y = 4,
             label = paste("nRMSE = ",round(nRMSE_8site[5], 
                                            digits = 2)), size = 4) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10), 
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# S. Fork Little Red River, AR
(gpp_plot8.4 <- ggplot(df_pred8$nwis_07075250, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Small Order - Turbulent Flow - Low Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2016-06-01"), y = 2.75,
             label = paste("nRMSE = ",round(nRMSE_8site[6], 
                                            digits = 2)), size = 4) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10), 
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# Red River, ND
(gpp_plot8.5 <- ggplot(df_pred8$nwis_05082500, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Large Order - Steady Flow - High Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2016-10-15"), y = 9,
             label = paste("nRMSE = ",round(nRMSE_8site[4], 
                                            digits = 2)), size = 4) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10), 
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# Snake River, WY
(gpp_plot8.6 <- ggplot(df_pred8$nwis_13013650, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median, group = seq), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Large Order - Steady Flow - Low Light") +
    scale_x_date(date_labels = "%b %Y") +
    #scale_x_break(c(as.Date("2011-01-01"), as.Date("2014-01-01"))) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, group = seq),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2015-01-01"), y = 27,
             label = paste("nRMSE = ",round(nRMSE_8site[8], 
                                            digits = 2)), size = 4) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10), 
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# Raisin River, MI 
(gpp_plot8.7 <- ggplot(df_pred8$nwis_04176500, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Large Order - Turbulent Flow - High Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2014-03-01"), y = 45,
             label = paste("nRMSE = ",round(nRMSE_8site[3], 
                                            digits = 2)), size = 4) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10), 
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# Rio Grande, TX
(gpp_plot8.8 <- ggplot(df_pred8$nwis_08374550, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median, group = seq), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Large Order - Turbulent Flow - Low Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, group = seq),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2012-03-01"), y = 6,
             label = paste("nRMSE = ",round(nRMSE_8site[7], 
                                            digits = 2)), size = 4) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 10),
          axis.title.x = element_text(size=10), 
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10)))

# Combine and export.
(fig_nRMSE <- gpp_plot8.1 + gpp_plot8.2 + 
    gpp_plot8.3 + gpp_plot8.4 +
    gpp_plot8.5 + gpp_plot8.6 +
    gpp_plot8.7 + gpp_plot8.8 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 4))

# ggsave(fig_nRMSE,
#        filename = "figures/beartooth_spring23/nRMSE_8panel_051223.jpg",
#        width = 22,
#        height = 22,
#        units = "cm") # empty dates fixed :)

# End of script.
