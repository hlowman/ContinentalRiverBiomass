## Resilience of Stream Productivity to Disturbance
## October 12, 2022
## Heili Lowman

# This file will be used for the creation of additional plots and figures not
# part of the formal analysis workflow.

#### Setup #####

# Load packages.
lapply(c("tidyverse", "cowplot", "gt",
         "lubridate", "reshape2", "webshot2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", "here"), require, character.only=T)

# Load data required for figures below.

# Site information
site_info <- readRDS("data_working/NWIS_198sitesinfo_101222.rds")

# Unedited site information from hypoxia project.
hypox_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")

# Raw data used to fit models on Teton.
site_dat <- readRDS("data_working/df_182sites_Qmaxnorm_allSL.rds")

#### Summary Stats Table ####

# Join together site info and raw data.
site_all <- left_join(site_dat, site_info)

# Also join with all site info courtesy of hypoxia dataset.
site_all <- left_join(site_all, hypox_info, by = c("site_name" = "SiteID"))

# Going to make a gt summary table of the sites used for the run on 10/12/22.
# These will calculate mean values at the timescales denoted below.
site_summary <- site_all %>%
  group_by(site_name) %>%
  summarize(GPP_mean = mean(GPP, na.rm = TRUE), # daily
            temp_mean = mean(temp, na.rm = TRUE), # daily
            Q_mean = mean(Q, na.rm = TRUE), # daily
            light_mean = mean(PAR_surface, na.rm = TRUE), # daily
            #ws_mean = mean(NHD_AREASQKM, na.rm = TRUE), # by site
            records = n()) %>% # total
  ungroup()

# Pivot for easier creation of the gt table.
summary_pivot <- site_summary %>%
  pivot_longer(!site_name, names_to = "covariate", values_to = "value") %>%
  mutate(cov_f = factor(covariate, levels = c("records",
                                              "GPP_mean",
                                              "Q_mean",
                                              "light_mean",
                                              "temp_mean"))) %>%
  group_by(cov_f) %>%
  summarize(Minimum = round(min(value, na.rm = TRUE), digits = 2),
            Median = round(median(value, na.rm = TRUE), digits = 2), 
            Mean = round(mean(value, na.rm = TRUE), digits = 2), 
            Maximum = round(max(value, na.rm = TRUE), digits = 2)) %>%
  ungroup() %>%
  # Adding new column for easier formatting
  # With unicode versions of super/subscripts pasted in
  mutate(Covariate = c("Total Number of Daily Records",
                       "Mean Daily Gross Primary Production (gO₂ m⁻² d⁻¹)",
                       "Mean Daily Discharge (m³ s⁻¹)",
                       "Mean Daily Light Availability (µmol m⁻² s⁻¹)",
                       "Mean Daily Temperature (°C)"))

# Turns out the hypoxia dataset doesn't have the Miss. R. watershed area size,
# which was obviously our maximum but wasn't appearing as such, so I'm leaving
# watershed size off of the summary statistics table for now.

( summary_table <- summary_pivot %>%
    select(!cov_f) %>% # remove old covariate column
  # Base table creation/call
  gt(rowname_col = "Covariate") %>%
    # Add helpful (sub)titles
    tab_header(title = "Site-Level Stream Data Summary") )

# I tried doing cumulative annual means for light and GPP, but due to the 
# sometimes large differences in length of record, I decided to keep 
# everything at a daily scale to make things more comparable.

# Export table.
# gtsave(summary_table,
#        "summary_data_table_101222.png",
#        path = "figures/teton_fall22") 

#### Timeseries Length Appendix Figure ####

# Making a histogram of the lengths of the timeseries for an appendix figure.

site_lengths <- site_all %>%
  group_by(site_name) %>%
  summarize(days = n()) %>% # count the number of rows/days
  ungroup() %>%
  mutate(years = days/365)

ggplot(site_lengths, aes(x = years)) +
  geom_histogram() +
  theme_bw()

# Adding bins for more intuitive plotting.
site_lengths <- site_lengths %>%
  mutate(bin = factor(case_when(days < 180 ~ "3 to 6 months",
                         days >= 180 & days < 365 ~ "6 months to 1 year",
                         days >= 365 & days < 730 ~ "1 to 2 years",
                         days >= 730 & days < 1460 ~ "2 to 4 years",
                         days >= 1460 & days < 2190 ~ "4 to 6 years",
                         days >= 2190 & days < 2920 ~ "6 to 8 years",
                         days >= 2920 ~ "8 to 10 years"),
                      levels = c("3 to 6 months", "6 months to 1 year",
                                 "1 to 2 years", "2 to 4 years",
                                 "4 to 6 years", "6 to 8 years",
                                 "8 to 10 years")))

(ts_hist <- ggplot(site_lengths, aes(x = bin)) +
  geom_histogram(stat = "count") +
  labs(x = "Length of Timeseries",
       y = "Number of Sites") +
  theme_bw())

# And export for use in manuscript Google doc.
# ggsave(ts_hist,
#        filename = "figures/teton_fall22/TS_length_hist_020723.jpg",
#        width = 20,
#        height = 10,
#        units = "cm")

# And calculate median timeseries length for inclusion in methods.
median(site_lengths$days) # 873.5
            
# End of script.
