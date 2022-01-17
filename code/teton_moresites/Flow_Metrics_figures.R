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
# (6) amplitude - not included in this code
# (7) phase of the seasonal signal - not included in this code

# I will be using the dataset containing 207 sites generated in the script 
# code/teton_moresites/NWIS_RiverSelection.R.

## ADDENDUM: Following investigating the above metrics, I will also be examining
# other metrics to better delineate storm events, per a discussion with Joanna in
# Jan 2022. These investigations will be in the section titled "Additional Metrics".

#### Load-in ####

# Load packages
lapply(c("ggplot2","cowplot","lubridate",
         "data.table","patchwork", "here",
         "calecopal", "sf", "viridis", "mapproj", "moments",
         "corrplot", "gt", "webshot", "bstfun", 
         #  "EcoHydRology", - removed bc it causes problems with the tidyverse/ggplot
         "tidyverse"), require, character.only=T)

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

# Calculate flow statistics by site.
site_summary <- site_scaled %>%
  group_by(site_name) %>% # group by site
  summarize(maxQ = max(Q, na.rm = TRUE), # adding maximum Q for table in figure below
            meanQ = mean(Q, na.rm = TRUE), # (1) mean
            cvQ = (sd(Q, na.rm = TRUE)/mean(Q, na.rm = TRUE)), # (2) coefficient of variation
            skewQ = skewness(Q, na.rm = TRUE), # (3) skewness - positive = pulled right
            kurtQ = kurtosis(Q, na.rm = TRUE), # (4) kurtosis - positive = pointy/leptokurtic
            ar1Q = acf_print(scaleQ)) %>% # (5) AR(1) correlation coefficient
          # AQ = sqrt((a^2)+(b^2)) # (6) amplitude of seasonal signal - not included in this code
          # phiQ = atan(-a/b) # (7) phase shift of seasonal signal - not included in this code
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
# In plain English,
# when AR(1) = 1, flow is more persistent from one day to the next at a given site
# when AR(1) < 0.5, flow is variable from one day to the next at a given site
# Following a discussion with Joanna on 1/11/22, 0.8 appears to be a good cutoff
# at which to delineate sites that don't encounter changes in flow and might be 
# good sites at which to remove the P (persistence) term.

(fig_mean_ar1Q <- site_summary %>%
    mutate(log_meanQ = log10(meanQ)) %>%
    ggplot(aes(ar1Q, log_meanQ)) +
    geom_point() +
    theme_bw() +
    labs(x = "AR(1) Correlation Coefficient of Discharge (Q)",
         y = "Log of Mean Discharge (Q)")) 
# as mean Q increases, AR(1) cc increases 
# (more flow = more persistent)

(fig_cv_ar1Q <- site_summary %>%
    ggplot(aes(ar1Q, cvQ)) +
    geom_point() +
    theme_bw() +
    labs(x = "AR(1) Correlation Coefficient of Discharge (Q)",
         y = "Coefficient of Variation of Discharge (Q)")) 
# as cvQ increases, AR(1) cc decreases 
# (less flashy flow [i.e., lower cvQ] = more persistent day-to-day discharge)

# join with stream info data to examine if stream order has anything to do with it
site_summary_info <- join(site_summary, sitesjoin, by = "site_name")

(fig_order_ar1Q <- site_summary_info %>%
    ggplot(aes(ar1Q, NHD_STREAMORDE)) +
    geom_point(size = 3) +
    theme_bw() +
    scale_y_continuous(breaks = seq(0,10,2)) +
    labs(x = "AR(1) Correlation Coefficient of Discharge (Q)",
         y = "Stream Order")) 
# higher order streams have greater AR(1) a.k.a. more persistent discharge
# (more flow = more persistent)
# but it seems there are a few lower order streams that have high persistence
# suggesting all the more that AR(1) might be a good measure of flow, regardless
# stream size

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

# Following meeting with Joanna on 1/11/2022, creating another figure
# hopefully to be included as a supplementary figure in the manuscript
# to justify AR(1)-based decision-making. I'll be highlighting three
# sites in particular and providing their metrics in a table below.

site_summary_log <- site_summary %>%
  mutate(log_meanQ = log10(meanQ))

site_name3 <- c("nwis_05579630", "nwis_04199500", "nwis_01673000") # sites of interest in table

site_summary_log_label <- site_summary_log %>%
  mutate(site_labels = case_when(site_name %in% site_name3 ~ site_name,
                                 TRUE ~ ""))

(fig_supp1A_AR1 <- site_summary_log_label %>%
    ggplot(aes(x = ar1Q, y = cvQ)) +
    # highlight certain data points
    annotate(geom = "point", x = 0.8220972, y = 1.317274, color = "orange", size = 5) + # kickapoo
    annotate(geom = "point", x = 0.5112334, y = 2.269440, color = "orange", size = 5) + # vermilion
    annotate(geom = "point", x = 0.2438063, y = 5.120188, color = "orange", size = 5) + # pamunkey
    geom_point(aes(color = log_meanQ), size = 3) +
    scale_color_viridis() +
    # label certain data points
    geom_label(data = site_summary_log_label %>% filter(site_name %in% site_name3),
               aes(ar1Q, cvQ, label = site_name), size = 4,
                    hjust = 0.5, vjust = -0.5, fill = NA) +
    # note, geom_text_repel here does NOT play nice, particularly with the legends
    theme_bw() +
    labs(x = "AR(1) Correlation Coefficient of Discharge (Q)",
         y = "Coefficient of Variation of Discharge (Q)",
         color = "Log of Mean\nDischarge (Q)") +
    theme(legend.position = c(0.8, 0.75),  # move legend inside plot
          legend.background = element_rect(fill = NA,
                                           size = 0.5, linetype = "solid",
                                           color = "black"))) # edit legend aesthetics

(fig_supp1B_AR1 <- site_summary_log %>%
  ggplot(aes(x = ar1Q)) +
  geom_bar(color = "black", fill = "gray50") +
  scale_x_binned() +
  #geom_density(alpha = 0.2, fill = "gray48") +
  theme_bw() +
  labs(x = "AR(1) Correlation Coefficient of Discharge (Q)",
       y = "Site Count"))

# create data table to use in gt() below
site_summary_3 <- site_summary_log %>%
  filter(site_name %in% site_name3) %>%
  mutate(long_name = c("Kickapoo Creek, IL", "Vermilion River, OH", "Pamunkey River, VA")) %>%
  select(site_name, long_name, maxQ, meanQ, log_meanQ, cvQ, ar1Q)

(fig_supp1C_AR1 <- site_summary_3 %>%
    gt() %>% # create gt table of dataset above
    cols_label(
      site_name = "Site ID",
      long_name = "Site Name",
      maxQ = "Maximum",
      meanQ = "Mean",
      log_meanQ = "Log Mean",
      cvQ = "CV",
      ar1Q = "AR(1)") %>% # rename columns
    fmt_number(
      columns = 3:7,
      decimals = 2)) # rounds everything to 2 decimal places

# convert gt object into ggplot object using as_ggplot() function in bstfun package
fig_supp1C_AR1 <- as_ggplot(fig_supp1C_AR1)

# combine into a single plot
fig_supp1 <- fig_supp1A_AR1 + fig_supp1B_AR1 + fig_supp1C_AR1 + plot_layout(ncol = 2)

fig_supp1 + plot_annotation(tag_levels = "A")

# ggsave(("figures/teton_moresites/supp_fig_AR1.png"),
#        width = 24,
#        height = 16,
#        units = "cm"
# )

# Three storm comparison - trying to develop a conceptual figure to explain AR(1)
# reasoning with more direct relation to flow/discharge.

# Using the same three sites highlighted above (since the trends work so nicely)

# No storms - low AR(1): Kickapoo Creek, nwis_05579630
# Small storms - medium AR(1): Vermilion River, nwis_04199500
# Large storms - high AR(1): Pamunkey River, nwis_01673000

lowAR_site_dat <- site_subset %>%
  filter(site_name == "nwis_05579630") %>%
  mutate(month = month(date)) %>%
  filter(year == 2012 & month == 5)

# calculate AR(1) for this subset
# Using site_deseason and acf_print function created in the section above
lowAR_site <- site_scaled %>%
  filter(site_name == "nwis_05579630") %>%
  filter(year == 2012 & month == 5) %>%
  summarize(ar1Q = acf_print(scaleQ)) # 0.63

medAR_site_dat <- site_subset %>%
  filter(site_name == "nwis_04199500") %>%
  mutate(month = month(date)) %>%
  filter(year == 2014 & month == 5)

# calculate AR(1) for this subset
# Using site_deseason and acf_print function created in the section above
medAR_site <- site_scaled %>%
  filter(site_name == "nwis_04199500") %>%
  filter(year == 2014 & month == 5) %>%
  summarize(ar1Q = acf_print(scaleQ)) # 0.19

highAR_site_dat <- site_subset %>%
  filter(site_name == "nwis_01673000") %>%
  mutate(month = month(date)) %>%
  filter(year == 2014 & month == 5)

# calculate AR(1) for this subset
# Using site_deseason and acf_print function created in the section above
highAR_site <- site_scaled %>%
  filter(site_name == "nwis_01673000") %>%
  filter(year == 2014 & month == 5) %>%
  summarize(ar1Q = acf_print(scaleQ)) # 0.66

# Ok, so I think this may not be the best comparison (i.e., apples to apples)
# since this is across sites with varying monthly mean flow, so I'm going
# to zero in on the Pamunkey site and pick out no, low, and high storm events.

# rough plot to see which months to parse out
(pamunkey_2014 <- site_scaled %>%
  filter(site_name == "nwis_01673000") %>%
  filter(year == 2014) %>%
  ggplot() +
  geom_line(aes(x = date, y = Q)))

# again, using "site_scaled" dataset from above in which data has already been
# detrended against long-term monthly means and z-scored prior to applying the
# acf_print function

# february - small storms
pamunkey_feb_14_ar1 <- site_scaled %>%
  filter(site_name == "nwis_01673000") %>%
  filter(year == 2014 & month == 2) %>%
  summarize(ar1Q = acf_print(scaleQ)) # 0.67

# may - HUGE storms
pamunkey_may_14_ar1 <- site_scaled %>%
  filter(site_name == "nwis_01673000") %>%
  filter(year == 2014 & month == 5) %>%
  summarize(ar1Q = acf_print(scaleQ)) # 0.66

# august - no storms
pamunkey_aug_14_ar1 <- site_scaled %>%
  filter(site_name == "nwis_01673000") %>%
  filter(year == 2014 & month == 8) %>%
  summarize(ar1Q = acf_print(scaleQ)) # 0.50

# overall
pamunkey_14_ar1 <- site_scaled %>%
  filter(site_name == "nwis_01673000") %>%
  filter(year == 2014) %>%
  summarize(ar1Q = acf_print(scaleQ)) # 0.76

# Hmmm, maybe the AR(1) isn't the best linked with particular storms
(pamunkey_2014_ar1 <- site_scaled %>%
    filter(site_name == "nwis_01673000") %>%
    filter(year == 2014) %>%
    #mutate(Date = ymd(date)) %>%
    ggplot() +
    geom_line(aes(x = date, y = Q)) +
    annotate(xmin = ymd('2014-02-01'), xmax = ymd('2014-02-28'), ymin = -5, ymax = 500, 
             geom = "rect", fill = NA, color = "#69B9FA", size = 2) +
    annotate(xmin = ymd('2014-05-01'), xmax = ymd('2014-05-31'), ymin = -5, ymax = 500, 
             geom = "rect", fill = NA, color = "#4B8FF7", size = 2) +
    annotate(xmin = ymd('2014-08-01'), xmax = ymd('2014-08-31'), ymin = -5, ymax = 500, 
             geom = "rect", fill = NA, color = "#6B6D9F", size = 2) +
    annotate("text", x = ymd('2014-02-15'), y = 515, label = "AR(1) = 0.67") +
    annotate("text", x = ymd('2014-05-15'), y = 515, label = "AR(1) = 0.66") +
    annotate("text", x = ymd('2014-08-15'), y = 515, label = "AR(1) = 0.50") +
    annotate("text", x = ymd('2014-11-15'), y = 475, label = "2014\nAR(1) = 0.76") +
    annotate("text", x = ymd('2014-11-15'), y = 375, label = "Overall\nAR(1) = 0.82") +
    theme_bw() +
    labs(x = "Date",
         y = "Discharge (cm/s)"))

# This too isn't a great example. So switching back to trying to find three site-years
# with distinct storm regimes.

# calculate AR1 at site with few standout storms
lowcvsite_ar1 <- site_scaled %>%
  filter(site_name == "nwis_04137005") %>%
  filter(year == 2010) %>%
  summarize(ar1Q = acf_print(scaleQ)) # 0.92

# separate out data at site with few standout storms
lowcvsite_dat <- site_scaled %>%
  filter(site_name == "nwis_04137005") %>%
  filter(year == 2010) %>%
  mutate(ar1 = 0.92)

# calculate AR1 at site with a couple of storms
medcvsite_ar1 <- site_scaled %>%
  filter(site_name == "nwis_03067510") %>%
  filter(year == 2012) %>%
  summarize(ar1Q = acf_print(scaleQ)) # 0.52

# separate out data at site with a couple of storms
medcvsite_dat <- site_scaled %>%
  filter(site_name == "nwis_03067510") %>%
  filter(year == 2012) %>%
  mutate(ar1 = 0.52)

# calculate AR1 at site with one large storm
highcvsite_ar1 <- site_scaled %>%
  filter(site_name == "nwis_02336728") %>%
  filter(year == 2016) %>%
  summarize(ar1Q = acf_print(scaleQ)) # 0.04

# separate out data at site with one large storm
highcvsite_dat <- site_scaled %>%
  filter(site_name == "nwis_02336728") %>%
  filter(year == 2016) %>%
  mutate(ar1 = 0.04)

# join 3 sites of data together
ar1_3sites_bound <- bind_rows(lowcvsite_dat, medcvsite_dat, highcvsite_dat)

# create figure of data above

# create list for facet labels
f_labels <- c("0.04" = "Utoy Creek, GA - 2016\nAR(1) = 0.04",
  "0.52" = "Shavers Fork, WV - 2012\nAR(1) = 0.52",
  "0.92" = "Au Sable River,MI - 2010\nAR(1) = 0.92")

(fig_supp1_2 <- ar1_3sites_bound %>%
    mutate(ar1_f = factor(ar1)) %>%
    ggplot() +
    geom_line(aes(x = date, y = Q, color = ar1_f), size = 1.5) +
    scale_color_manual(values = c("#69B9FA", "#4B8FF7", "#6B6D9F")) +
    scale_x_date(date_labels = "%b") +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_rect(fill = NA)) +
    labs(x = "Date",
         y = "Discharge (cm/s)") +
    facet_wrap(.~ as.character(ar1_f), scales = "free", labeller = as_labeller(f_labels)))

# ggsave(("figures/teton_moresites/supp_fig_AR1_v2.png"),
#        width = 20,
#        height = 8,
#        units = "cm"
# )

#### Additional Metrics ####

# Just to make life a bit easier at the start, I'm going to start with just a few sites, some
# of which will match the 3 sites I've examined above.

# [nwis_02336728] Utoy Creek (GA) 2016 - one larger storm
# [nwis_03067510] Shavers Fork (WV) 2012 - a few storms
# [nwis_04137005] Au Sable River (MI) 2010 - multiple smaller storms
# [nwis_05579630] Kickapoo Creek (IL) 2012 - no storms

# also:
# [nwis_04085108] East River (WI) 2016 - no storms (summer only)

list_of_4 <- c("nwis_02336728", "nwis_03067510", "nwis_04137005", "nwis_05579630")

site_subset_4 <- site_subset %>%
  filter(site_name %in% list_of_4)

# First, I'll examine if the range > mean. If false, we would drop the P term
utoy_creek_m1 <- site_subset %>%
  filter(site_name == "nwis_02336728" & year == 2016) %>%
  summarize(rangeQ = max(range(Q, na.rm = TRUE)) - min(range(Q, na.rm = TRUE)),
            meanQ = mean(Q, na.rm = TRUE)) # Range > Mean == TRUE

shavers_fork_m1 <- site_subset_4 %>%
  filter(site_name == "nwis_03067510" & year == 2012) %>%
  summarize(rangeQ = max(range(Q, na.rm = TRUE)) - min(range(Q, na.rm = TRUE)),
            meanQ = mean(Q, na.rm = TRUE)) # Range > Mean == TRUE

au_sable_m1 <- site_subset_4 %>%
  filter(site_name == "nwis_04137005" & year == 2010) %>%
  summarize(rangeQ = max(range(Q, na.rm = TRUE)) - min(range(Q, na.rm = TRUE)),
            meanQ = mean(Q, na.rm = TRUE)) # Range > Mean == TRUE (but it's close)

kickapoo_creek_m1 <- site_subset_4 %>%
  filter(site_name == "nwis_05579630" & year == 2012) %>%
  summarize(rangeQ = max(range(Q, na.rm = TRUE)) - min(range(Q, na.rm = TRUE)),
            meanQ = mean(Q, na.rm = TRUE)) # Range > Mean == TRUE (but it's VERY CLOSE)

east_river_m1 <- site_subset %>%
  filter(site_name == "nwis_04085108" & year == 2016 & doy >=153 & doy <= 244) %>%
  summarize(rangeQ = max(range(Q, na.rm = TRUE)) - min(range(Q, na.rm = TRUE)),
            meanQ = mean(Q, na.rm = TRUE)) # Range > Mean == TRUE (but it's SUPER DUPER close)

# datasets for easier use

utoy_dat <- site_subset_4 %>%
  filter(site_name == "nwis_02336728" & year == 2016)

shavers_dat <- site_subset_4 %>%
  filter(site_name == "nwis_03067510" & year == 2012)

au_sable_dat <- site_subset_4 %>%
  filter(site_name == "nwis_04137005" & year == 2010)

kickapoo_dat <- site_subset_4 %>%
  filter(site_name == "nwis_05579630" & year == 2012)

# east_dat <- site_subset_4 %>%
#   filter(site_name == "nwis_04085108" & year == 2016 & doy >=153 & doy <= 244)

# Next, I'm going to see if flow exceeds 150% of mean flow at any point
utoy_mean <- 0.81
utoy_creek_m2 <- utoy_dat %>%
  mutate(exceedance = case_when(Q > (1.5*utoy_mean) ~ 1,
                                Q <= (1.5*utoy_mean) ~ 0))
sum(utoy_creek_m2$exceedance) # 43 days

shavers_mean <- 3.41
shavers_m2 <- shavers_dat %>%
  mutate(exceedance = case_when(Q > (1.5*shavers_mean) ~ 1,
                                Q <= (1.5*shavers_mean) ~ 0))
sum(shavers_m2$exceedance) # 57 days

au_sable_mean <- 31.01
au_sable_m2 <- au_sable_dat %>%
  mutate(exceedance = case_when(Q > (1.5*au_sable_mean) ~ 1,
                                Q <= (1.5*au_sable_mean) ~ 0))
sum(au_sable_m2$exceedance) # 7 days

kickapoo_mean <- 0.08
kickapoo_m2 <- kickapoo_dat %>%
  mutate(exceedance = case_when(Q > (1.5*kickapoo_mean) ~ 1,
                                Q <= (1.5*kickapoo_mean) ~ 0))
sum(kickapoo_m2$exceedance) # 45 days

# Now, I'm going to investigate the baseflow function in the EcoHydRology package.

# get an approximation for baseflow using a 3 pass filter and a value of 0.925
# rturns a 2 column data frame with first column - baseflow and second column -
# quickflow, in same units as input

# library(EcoHydRology)
#
# utoy_bfs <- BaseflowSeparation(utoy_dat$Q, passes = 3)
# 
# utoy_together <- cbind(utoy_dat, utoy_bfs)
# 
# ggplot(utoy_together) +
#   geom_line(aes(x = date, y = Q), color = "blue") +
#   geom_line(aes(x = date, y = bt), colo = "black") +
#   labs(x = "Date", y = "Discharge (cm/s)")
# 
# shavers_bfs <- BaseflowSeparation(shavers_dat$Q, passes = 3)
# shavers_together <- cbind(shavers_dat, shavers_bfs)
# 
# au_sable_bfs <- BaseflowSeparation(au_sable_dat$Q, passes = 3)
# au_sable_together <- cbind(au_sable_dat, au_sable_bfs)
# 
# kickapoo_bfs <- BaseflowSeparation(kickapoo_dat$Q, passes = 3)
# kickapoo_together <- cbind(kickapoo_dat, kickapoo_bfs)

# dataset for plotting

# fm_4sites_bound <- bind_rows(utoy_together, 
                             # shavers_together, 
                             # au_sable_together, 
                             # kickapoo_together)

# write_rds(fm_4sites_bound, "data_working/flowmetrics_baseflow_4sites_01172022.rds")

fm_4sites_bound <- readRDS("data_working/flowmetrics_baseflow_4sites_01172022.rds")

# create figure of data above

# create list for facet labels
f_labels2 <- c("nwis_02336728" = "Utoy Creek, GA - 2016\nRange > Mean TRUE\nExceedance of 150% Mean 43 days",
              "nwis_03067510" = "Shavers Fork, WV - 2012\nRange > Mean TRUE\nExceedance of 150% Mean 57 days",
              "nwis_04137005" = "Au Sable River, MI - 2010\nRange > Mean TRUE\nExceedance of 150% Mean 7 days",
              "nwis_05579630" = "Kickapoo Creek, IL - 2012\nRange > Mean TRUE\nExceedance of 150% Mean 45 days")

(fig_supp1_3 <- fm_4sites_bound %>%
    ggplot() +
    geom_line(aes(x = date, y = Q, color = site_name), size = 1.5) +
    geom_line(aes(x = date, y = bt), color = "black", size = 1) + # baseflow in black
    scale_color_manual(values = c("#69B9FA", "#4B8FF7", "#6B6D9F", "#D46F10")) +
    scale_x_date(date_labels = "%b") +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_rect(fill = NA)) +
    labs(x = "Date",
         y = "Discharge (cm/s)") +
    facet_wrap(.~ as.character(site_name), scales = "free", labeller = as_labeller(f_labels2)))

# ggsave(("figures/teton_moresites/supp_fig_AR1_v3.png"),
#        width = 20,
#        height = 20,
#        units = "cm"
# )

# End of script.
