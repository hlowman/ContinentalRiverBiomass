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

# Raw data used to fit models on Teton.
site_dat <- readRDS("data_working/df_182sites_Qmaxnorm_allSL.rds")

#### Summary Stats Table ####

# Join together site info and raw data.
site_all <- left_join(site_dat, site_info)

# Going to make a gt summary table of the sites used for the run on 10/12/22.
# These will calculate mean values at the timescales denoted below.
site_summary <- site_all %>%
  group_by(site_name) %>%
  summarize(GPP_mean = mean(GPP, na.rm = TRUE), # daily
            temp_mean = mean(temp, na.rm = TRUE), # daily
            Q_mean = mean(Q, na.rm = TRUE), # daily
            light_mean = mean(PAR_surface, na.rm = TRUE), # daily
            records = n()) %>% # total
  ungroup()

# Pivot for easier creation of the gt table.
summary_pivot <- site_summary %>%
  pivot_longer(!site_name, names_to = "covariate", values_to = "value") %>%
  group_by(covariate) %>%
  summarize(Minimum = round(min(value, na.rm = TRUE), digits = 2),
            Median = round(median(value, na.rm = TRUE), digits = 2), 
            Mean = round(mean(value, na.rm = TRUE), digits = 2), 
            Maximum = round(max(value, na.rm = TRUE), digits = 2)) %>%
  ungroup() %>%
  # Adding new column for easier formatting
  # With unicode versions of super/subscripts pasted in
  mutate(Covariate = c("Mean Daily Gross Primary Production (gO₂ m⁻² d⁻¹)",
                       "Mean Daily Light Availability (µmol m⁻² s⁻¹)",
                       "Mean Daily Discharge (m³ s⁻¹)",
                       "Total Number of Daily Records",
                       "Mean Daily Temperature (°C)"))

( summary_table <- summary_pivot %>%
    select(!covariate) %>% # remove old covariate column
  # Base table creation/call
  gt(rowname_col = "Covariate") %>%
    # Add helpful (sub)titles
    tab_header(title = "Site-Level Stream Data Summary") )

# I tried doing cumulative annual means for light and GPP, but due to the 
# sometimes large differences in length of record, I decided to keep 
# everything at a daily scale to make things more comparable.

# Export table.
gtsave(summary_table,
       "summary_data_table_101222.png",
       path = "figures/teton_fall22") 
            
# End of script.
