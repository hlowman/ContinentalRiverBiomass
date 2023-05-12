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
ggsave(fig_yield_med,
       filename = "figures/beartooth_spring23/amax_11panel_051223.jpg",
       width = 22,
       height = 26,
       units = "cm") # n = 152

#### Qc Appendix Figure ####

# Distribution of c values:
(fig0c <- ggplot(dat_out_full_141, aes(x = c_med)) +
   geom_histogram(bins = 60, alpha = 0.8, 
                  fill = "#262E43", color = "#262E43") +
   labs(x = expression(Critical~Disturbance~Threshold~(Q[c])),
        y = "Count") +
   theme_bw())

# Examining uncertainty
(fig0c2 <- ggplot(dat_out_full_141, aes(x = c_med, y = site_name)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#405F8A") +
    geom_linerange(alpha = 0.8, 
                   color = "#405F8A",
                   aes(xmin = `2.5%`, xmax = `97.5%`)) +
    #scale_x_log10() +
    #scale_y_log10() + 
    labs(x = expression(c)) +
    theme_bw() +
    theme(axis.text.y=element_blank()))

# Add a new column to quantify range of uncertainty
dat_out_full_141 <- dat_out_full_141 %>%
  mutate(range_c = `97.5%` - `2.5%`)

plot(dat_out_full_141$c_med, dat_out_full_141$range_c)

# Distribution of Qc/Q2 values:
(fig1qcq2 <- ggplot(dat_out_full_141, aes(x = Qc_Q2yr)) +
    geom_histogram(bins = 60, alpha = 0.8, 
                   fill = "#262E43", color = "#262E43") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    labs(x = expression(Q[c]:Q[2~yr]),
         y = "Count") +
    theme_bw())

# CV of Discharge vs. Qc:Q2: note, x axis on LOG SCALE
(fig2qcq2 <- ggplot(dat_out_full_141, aes(x = cvQ, y = Qc_Q2yr)) +
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

(fig2c2 <- ggplot(dat_out_full_141, aes(x = cvQ, y = c_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#405F8A") +
    geom_linerange(alpha = 0.8, 
                   color = "#405F8A",
                   aes(ymin = `2.5%`, ymax = `97.5%`)) +
    #scale_x_log10() +
    #scale_y_log10() + 
    labs(x = expression(CV[Q]),
         y = expression(c)) +
    theme_bw())

# Mean Daily Light Availability vs. c:
(fig3qcq2 <- ggplot(dat_out_full_141, aes(x = meanL, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.8, size = 3,
               color = "#E6A45A") +
    labs(x = expression(Mean~Daily~PAR),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Stream Order vs. c: Removing singular site w/o order info for now.
(fig4qcq2 <- ggplot(dat_out_full_141 %>%
                      filter(!is.na(Order)), aes(x = Order, y = Qc_Q2yr)) +
    geom_boxplot(alpha = 0.6, color = "black", fill = "#5792CC") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() +
    labs(x = expression(Stream~Order),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Stream Width vs. c: note, x axis LOG SCALED
(fig4.1qcq2 <- ggplot(dat_out_full_141, aes(x = width_med, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#1E2F46") +
    geom_linerange(alpha = 0.8, 
                   color = "#1E2F46",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_log10() + 
    scale_y_log10() + 
    labs(x = expression(River~Width~(m)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Latitude vs. c:
(fig5qcq2 <- ggplot(dat_out_full_141, aes(x = Lat_WGS84, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#7E8C69") +
    geom_linerange(alpha = 0.8, 
                   color = "#7E8C69",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(x = expression(Latitude),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Longitude vs. c:
(fig6qcq2 <- ggplot(dat_out_full_141, aes(x = Lon_WGS84, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.6, size = 3, color = "#E38678") +
    geom_linerange(alpha = 0.8, 
                   color = "#E38678",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    labs(x = expression(Longitude),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Additional exploratory plot:
(figextra <- ggplot(dat_out_full_141, aes(x = Lon_WGS84, y = cvQ)) +
    geom_point(alpha = 0.6, size = 3, color = "black") +
    labs(x = expression(Longitude),
         y = expression(CV[Q])) +
    theme_bw())

# Catchment size vs. c: note, missing Miss. R.
# x-axis also LOG SCALED
(fig7qcq2 <- ggplot(dat_out_full_141, aes(x = NHD_AREASQKM, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#3B7D6E") +
    geom_linerange(alpha = 0.8, 
                   color = "#3B7D6E",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Watershed~Area~(km^2)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Land-use vs. c:
(fig8qcq2 <- ggplot(dat_out_full_141, aes(x = LU_category, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_boxplot(alpha = 0.6, color = "#6D4847", fill = "#6D4847") +
    labs(x = expression(Land~Use),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Mean daily GPP vs. Qc:Q2: X axis LOG SCALED
(fig9qcq2 <- ggplot(dat_out_full_141, aes(x = meanGPP, y = Qc_Q2yr)) +
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

# More Land Use vs. c:
(fig10qcq2 <- ggplot(dat_out_full_141, aes(x = NHD_RdDensCat, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#9EB45F") +
    geom_linerange(alpha = 0.8, 
                   color = "#9EB45F",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() + 
    labs(x = expression(Road~Density~by~Catchment~(km/km^2)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

(fig11qcq2 <- ggplot(dat_out_full_141, aes(x = NHD_RdDensWs, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#304969") +
    geom_linerange(alpha = 0.8, 
                   color = "#304969",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() + 
    labs(x = expression(Road~Density~by~Watershed~(km/km^2)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

(fig12qcq2 <- ggplot(dat_out_full_141, 
                     aes(x = NHD_PctImp2011Cat, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#E3907B") +
    geom_linerange(alpha = 0.8, 
                   color = "#E3907B",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(x = expression(Percent~Impervious~by~Catchment),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

(fig13qcq2 <- ggplot(dat_out_full_141, 
                     aes(x = NHD_PctImp2011Ws, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.6, size = 3, color = "#6D4847") +
    labs(x = expression(Percent~Impervious~by~Watershed),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Effect of Dams
# "95 indicates the least probable interference from a structure of a given type"
(fig14qcq2 <- ggplot(dat_out_full_141 %>%
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

# Exploring uncertainty in Qc estimates by dam categories
(fig14qcq2.2 <- ggplot(dat_out_full_141, aes(x = Dam, y = Qc_Q2yr)) +
    geom_linerange(alpha = 0.8, position = position_jitter(width = 0.15),
                   color = "#4D5B90",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() + 
    labs(x = expression(Likelihood~of~Influence~by~Dams),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Effect of Canals
# "95 indicates the least probable interference from a structure of a given type"
(fig15qcq2 <- ggplot(dat_out_full_141, aes(x = Canal, y = Qc_Q2yr)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#E59D7F", color = "#E59D7F") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(x = expression(Likelihood~of~Influence~by~Canals),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Precip
df_precip_Qc_141 <- left_join(dat_out_full_141, site_precip, by = c("site_name" = "SiteID"))

(fig16qcq2 <- ggplot(df_precip_Qc_141, aes(x = pre_mm_cyr, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#4D5B75") +
    geom_linerange(alpha = 0.8, 
                   color = "#4D5B75",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() +
    labs(x = expression(Mean~Annual~Precipitation~by~Catchment~(mm)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# HUC vs. Qc:Q2: 
(fig17qcq2 <- ggplot(dat_out_full_141, aes(x = huc2_id, y = Qc_Q2yr)) +
    geom_boxplot(alpha = 0.6, color = "black", fill = "#38557A") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() +
    labs(x = expression(Regional~Hydrological~Unit~Code),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Combine figures above and export for supplemental figure.
(fig_qcq2_supp <- fig9qcq2 + fig2qcq2 + fig17qcq2 +
    fig11qcq2 + fig14qcq2 + fig4.1qcq2 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 2))

# ggsave(fig_qcq2_supp,
#        filename = "figures/teton_fall22/QcQ2_6panel_022823.jpg",
#        width = 30,
#        height = 20,
#        units = "cm") # n = 141

# End of script.
