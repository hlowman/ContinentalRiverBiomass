## Measures of flow calculations
## January 5, 2022
## Heili Lowman

# The following script is intended to calculate the measures of flow laid out in
# Archfield et al. 2014 in order to better inform the grouping of sites that will
# and will not include a persistence (P) term in their applied Ricker model.

# The flow metrics will be calculated according to the manuscript's instructions
# and they include: (1) mean, (2) coefficient of variation, (3) skewness, (4) kurtosis,
# (5) autoregressive lag-one correlation coefficient, (6) amplitude, and (7) phase of the
# seasonal signal.

# I will be using the dataset containing 207 sites generated in the script 
# code/teton_moresites/NWIS_RiverSelection.R.

# Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork", "here",
         "calecopal", "sf", "viridis", "mapproj"), require, character.only=T)

# Double check working directory
here()

# Load datasets.
# Main dataset containing site information from Appling + Koenig
sitesjoin <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")

# Dataset with data itself
site_subset <- readRDS("data_working/NWIS_207sites_subset.rds")

# Note: Discharge here is in cm/s.

# Calculate flow statistics by site.
site_summary <- site_subset %>%
  group_by(site_name) %>% # group by site
  summarize(meanQ = mean(Q, na.rm = TRUE), # (1) mean
            cvQ = (sd(Q, na.rm = TRUE)/mean(Q, na.rm = TRUE))) %>% # (2) coefficient of variation
  ungroup() # don't forget it!!

# Plot results.
# Mean discharge
(fig_meanQ <- site_summary %>%
    mutate(log_meanQ = log10(meanQ)) %>%
    ggplot(aes(log_meanQ, site_name)) +
    geom_point() +
    theme_bw() +
    labs(x = "Log of Mean Discharge (Q)",
       y = "Site"))

new <- site_summary %>%
    mutate(log_meanQ = log10(meanQ))

hist(new$log_meanQ)

# Coefficient of variation of discharge
(fig_cvQ <- site_summary %>%
    #mutate(log_meanQ = log10(meanQ)) %>%
    ggplot(aes(cvQ, site_name)) +
    geom_point() +
    theme_bw() +
    labs(x = "Coefficient of Variation of Discharge (Q)",
         y = "Site"))

hist(site_summary$cvQ)

# End of script.