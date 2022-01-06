## Measures of flow calculations
## January 5, 2022
## Heili Lowman

# The following script is intended to calculate the measures of flow laid out in
# Archfield et al. 2014 in order to better inform the grouping of sites that will
# and will not include a persistence (P) term in their applied Ricker model.

# The flow metrics will be calculated according to the manuscript's instructions
# and they include: 
# (1) mean
# (2) coefficient of variation
# (3) skewness
# (4) kurtosis
# (5) autoregressive lag-one correlation coefficient
# (6) amplitude
# (7) phase of the seasonal signal

# I will be using the dataset containing 207 sites generated in the script 
# code/teton_moresites/NWIS_RiverSelection.R.

#### Load-in ####

# Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork", "here",
         "calecopal", "sf", "viridis", "mapproj", "moments",
         "corrplot"), require, character.only=T)

# Double check working directory
here()

# Load datasets.
# Main dataset containing site information from Appling + Koenig
sitesjoin <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")

# Dataset with data itself
site_subset <- readRDS("data_working/NWIS_207sites_subset.rds")

# Note: Discharge here is in cm/s.

#### Calculate Flow Statistics ####

# Additional calculations required to determine autoregressive lag-one correlation coefficient.
# Need to first deseasonalize data by subtracting long-term monthly means from all values.

# calculate long-term monthly means
monthly_means <- site_subset %>%
  mutate(month = month(date)) %>%
  group_by(site_name, month) %>%
  summarize(lt_meanQ = mean(Q, na.RM = TRUE)) %>%
  ungroup()

# calculate deseasonalized discharge
site_month <- site_subset %>%
  mutate(month = month(date))

site_deseason <- site_month %>%
  join(monthly_means, by = c("site_name", "month")) %>%
  mutate(deseasonQ = Q - lt_meanQ)

# Then, standardize all values to have a zero mean and unit variance.
# Can use the scale() function here, but I'm going to do it by hand
# and calculate the z-score: (x - mean(x)) / sd(x)
site_scaled <- site_deseason %>%
  mutate(scaleQ = (deseasonQ - mean(deseasonQ, na.rm = TRUE))/sd(deseasonQ, na.rm = TRUE))

# create function to pull out AR(1) correlation coefficient
acf_print <- function(x) {
  a1 <- acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)
  a1$acf[2]
}

# test at a single site to make sure it works
site_test <- site_scaled %>%
  filter(site_name == "nwis_03081000") %>%
  summarize(ar1Q = acf_print(scaleQ)) # yay :)

# Insert calculations for qzt function here.

# Calculate flow statistics by site.
site_summary <- site_scaled %>%
  group_by(site_name) %>% # group by site
  summarize(meanQ = mean(Q, na.rm = TRUE), # (1) mean
            cvQ = (sd(Q, na.rm = TRUE)/mean(Q, na.rm = TRUE)), # (2) coefficient of variation
            skewQ = skewness(Q, na.rm = TRUE), # (3) skewness - positive = pulled right
            kurtQ = kurtosis(Q, na.rm = TRUE), # (4) kurtosis - positive = pointy/leptokurtic
            ar1Q = acf_print(scaleQ)) %>% # (5) AR(1) correlation coefficient
          # AQ = sqrt((a^2)+(b^2)) # (6) amplitude of seasonal signal
          # phiQ = atan(-a/b) # (7) phase shift of seasonal signal
  ungroup() # don't forget it!!

#### Figures ####

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
    ggplot(aes(cvQ, site_name)) +
    geom_point() +
    theme_bw() +
    labs(x = "Coefficient of Variation of Discharge (Q)",
         y = "Site"))

hist(site_summary$cvQ)

(fig_mean_cvQ <- site_summary %>%
    mutate(log_meanQ = log10(meanQ)) %>%
    ggplot(aes(cvQ, log_meanQ)) +
    geom_point() +
    theme_bw() +
    labs(x = "Coefficient of Variation of Discharge (Q)",
         y = "Log of Mean Discharge (Q)")) # as mean Q increases, cvQ decreases (i.e., less flashy)

# Skewness of discharge
(fig_skewQ <- site_summary %>%
    ggplot(aes(skewQ, site_name)) +
    geom_point() +
    theme_bw() +
    labs(x = "Skewness of Discharge (Q)",
         y = "Site"))

hist(site_summary$skewQ) # all right-skewed

(fig_mean_skewQ <- site_summary %>%
    mutate(log_meanQ = log10(meanQ)) %>%
    ggplot(aes(skewQ, log_meanQ)) +
    geom_point() +
    theme_bw() +
    labs(x = "Skewness of Discharge (Q)",
         y = "Log of Mean Discharge (Q)")) # as mean Q increases, right-skewness decreases

# Kurtosis of discharge
(fig_kurtQ <- site_summary %>%
    ggplot(aes(kurtQ, site_name)) +
    geom_point() +
    theme_bw() +
    labs(x = "Kurtosis of Discharge (Q)",
         y = "Site"))

hist(site_summary$skewQ) # all leptokurtic

(fig_mean_kurtQ <- site_summary %>%
    mutate(log_meanQ = log10(meanQ)) %>%
    ggplot(aes(kurtQ, log_meanQ)) +
    geom_point() +
    theme_bw() +
    labs(x = "Kurtosis of Discharge (Q)",
         y = "Log of Mean Discharge (Q)")) # as mean Q increases, kurtosis decreases

# AR(1) of discharge
(fig_ar1Q <- site_summary %>%
    ggplot(aes(ar1Q, site_name)) +
    geom_point() +
    theme_bw() +
    labs(x = "AR(1) Correlation Coefficient of Discharge (Q)",
         y = "Site"))

hist(site_summary$ar1Q) # streamflow relatively persistent from 1 day to the next

(fig_mean_ar1Q <- site_summary %>%
    mutate(log_meanQ = log10(meanQ)) %>%
    ggplot(aes(ar1Q, log_meanQ)) +
    geom_point() +
    theme_bw() +
    labs(x = "AR(1) Correlation Coefficient of Discharge (Q)",
         y = "Log of Mean Discharge (Q)")) # as mean Q increases, AR(1) cc increases (more flow = more persistent)

# full correlation figure
for_corr <- site_summary %>%
  select(-site_name)

full_corr <- cor(for_corr)

corrplot(full_corr, # uses corrplot package
         tl.srt = 45, # tilts labels 45 degrees
         bg = "White", # color scheme/white background 
         title = "Correlation Plot of Flow Metrics", # adds title
         addCoef.col = "black", # adds coefficients
         type = "lower") # only the lower half

# End of script.
