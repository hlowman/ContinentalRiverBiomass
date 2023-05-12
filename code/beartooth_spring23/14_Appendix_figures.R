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
         "reshape2","ggExtra"), 
       require, character.only=T)

#### Data ####

# Import necessary datasets.

# First, the data for the maximum accrual (amax) models.
dat_amax <- readRDS("data_working/amax_covariates_152sites_051123.rds")

# Next, the data for the Qc:Q2yr models.
dat_Qc <- readRDS("data_working/Qc_covariates_138sites_051123.rds")

# And, the original data used in modeling.
dat_in <- readRDS("data_working/df_181sites_Qmaxnorm_SavoySL.rds")

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

# End of script.
