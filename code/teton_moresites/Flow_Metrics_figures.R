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
         y = "Log of Mean Discharge (Q)")) # as mean Q increases, AR(1) cc increases 
# (more flow = more persistent)

(fig_cv_ar1Q <- site_summary %>%
    ggplot(aes(ar1Q, cvQ)) +
    geom_point() +
    theme_bw() +
    labs(x = "AR(1) Correlation Coefficient of Discharge (Q)",
         y = "Coefficient of Variation of Discharge (Q)")) # as cvQ increases, AR(1) cc decreases 
# (less flashy flow [i.e., lower cvQ] = more persistent day-to-day discharge)

# join with stream info data to examine if stream order has anything to do with it
site_summary_info <- join(site_summary, sitesjoin, by = "site_name")

(fig_order_ar1Q <- site_summary_info %>%
    ggplot(aes(ar1Q, NHD_STREAMORDE)) +
    geom_point(size = 3) +
    theme_bw() +
    scale_y_continuous(breaks = seq(0,10,2)) +
    labs(x = "AR(1) Correlation Coefficient of Discharge (Q)",
         y = "Stream Order")) # higher order streams have greater AR(1) a.k.a. more persistent discharge

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

# Based on persistence curves, choosing 3 sites where P doesn't vary and 3 sites where P responds
# at varying "steepnesses" to discharge to examine how these sites' flow metrics vary.

# Filter dataset for desired sites.
sites <- c("nwis_023362095", "nwis_03298250", "nwis_05406457", "nwis_01632900", "nwis_04125460", "nwis_04121944")

rsites <- c("nwis_01632900", "nwis_04125460", "nwis_04121944")

nrsites <- c("nwis_023362095", "nwis_03298250", "nwis_05406457")

site_select <- site_summary %>%
  filter(site_name %in% sites) %>%
  pivot_longer(!site_name, names_to = "metric", values_to = "value") %>%
  mutate(category = case_when(site_name %in% rsites ~ "responsive",
                              site_name %in% nrsites ~ "non-responsive"))

(fig_paneled <- site_select %>%
  ggplot(aes(value, site_name, color = category)) +
  geom_point(size = 3) +
  scale_color_manual(values = cal_palette("seagrass")) +
  theme_bw() +
  facet_wrap(.~metric, scales = "free"))

# Not super informative, so going to examine these values in the context of critical discharge instead.

# Import results from previous 34 site run.
run1_params <- readRDS("data_working/teton_34rivers_model_parameters_090821.rds")

site_withc <- run1_params %>%
  select(site_name, c_mean) %>% # chose only columns of interest
  join(site_summary, by = "site_name") %>% # join with table from above (dropping all but 34 sites)
  pivot_longer(cols = meanQ:ar1Q, names_to = "metric", values_to = "value")

(fig_paneled2 <- site_withc %>%
    ggplot(aes(value, site_name, color = c_mean)) +
    geom_point(size = 3) +
    scale_color_viridis() +
    theme_bw() +
    theme(axis.text.y = element_blank()) +
    facet_wrap(.~metric, scales = "free"))

# Another exploration based on mean Q.

site_select2 <- site_summary %>%
  mutate(log_meanQ = log10(meanQ)) %>%
  select(-meanQ) %>%
  pivot_longer(cols = cvQ:ar1Q, names_to = "metric", values_to = "value")

(fig_paneled3 <- site_select2 %>%
    ggplot(aes(value, site_name, color = log_meanQ)) +
    geom_point(size = 3) +
    scale_color_viridis() +
    theme_bw() +
    theme(axis.text.y = element_blank()) +
    facet_wrap(.~metric, scales = "free"))

# End of script.
