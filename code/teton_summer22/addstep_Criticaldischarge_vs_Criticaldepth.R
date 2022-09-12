## Resilience of Stream Productivity to Disturbance
## August 31, 2022
## Heili Lowman

# The following set of scripts will walk through the steps necessary to
# prep and send data to Teton as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_raw" and "data_working" folders have been ignored
# using git.ignore, so links to the raw data sets are provided in the step1
# file. If you are accessing the code via GitHub, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running
# the code.

# NOTE: This script is specifically designed to compared the Ricker model's
# calculations of "c" to critical discharge calculated based on
# geomorphic data provided by Jud Harvey.

#### Setup ####

## Load packages
lapply(c("tidyverse", "lubridate", "here", 
         "data.table", "pipeR", "patchwork",
         "ggplot2"), require, character.only=T)

#### Data Import ####

# Load site-level info (hypoxia and Appling datasets).
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
site <- fread("data_raw/site_data.tsv")

# Load site-level parameters dataset with all iterations included.
dat_out_q10 <- readRDS("data_working/teton_190rivers_model_site_params_all_iterations_082422.rds")
dat_out_max <- readRDS("data_working/teton_4rivers_model_site_params_all_iterations_090822.rds")
dat_out_max2 <- readRDS("data_working/teton_2rivers_model_site_params_all_iterations_091222.rds")

# Load in original data fed into the model.
dat_in <- readRDS("data_working/df_190sites_10yrQnorm_allSL.rds")
dat_in2 <- readRDS("data_working/list_4sites_meanmaxQnorm_allSL.rds")
dat_in2 <- map_df(dat_in2, ~as.data.frame(.x), .id="site_name")
dat_in3 <- readRDS("data_working/list_2sites_meanmaxQnorm_allSL.rds")
dat_in3 <- map_df(dat_in3, ~as.data.frame(.x), .id="site_name")

#### Data Selection ####

# Geomorphic data is only available for the following sites:
# nwis_05451210: South Fork Iowa River, New Providence, IA
# nwis_05524500: Iroquois River, Foresman, IN
# nwis_05579630: Kickapoo Creek, Bloomington, IL
# nwis_01645704: Difficult Run, Fairfax, VA
# nwis_0165389205: Accotink Creek, Ranger Road, VA
# nwis_05515500: Kankakee River, Davis, IN

# Filtering by 4 sites
my_list <- c("nwis_05451210","nwis_05524500",
             "nwis_05579630","nwis_01645704")

# Filtering by all 6 sites
my_list6 <- c("nwis_05451210","nwis_05524500",
              "nwis_05579630","nwis_01645704",
              "nwis_0165389205", "nwis_05515500")

dat_out_q10_6 <- dat_out_q10 %>%
  filter(site_name %in% my_list6)

# And calculate normalized c statistic values
c_summary_q10 <- dat_out_q10_6 %>%
  group_by(site_name) %>%
  summarize(mean_c = mean(c),
            med_c = median(c),
            q2.5_c = quantile(c, c(.025)),
            q97.5_c = quantile(c, c(.975))) %>%
  ungroup()

# Identify ten year flood values to back calculate c in cfs.
Q10_6 <- dat_in %>%
  filter(site_name %in% my_list6) %>%
  distinct(site_name,RI_10yr_Q,RI_10yr_Q_cms)

# Join these datasets.
c_Q10_6 <- left_join(c_summary_q10, Q10_6) # joined by "site_name"

# Calculate "c" in cubic meters per second.
c_Q10_6 <- c_Q10_6 %>%
  mutate(mean_c_cms = mean_c * RI_10yr_Q_cms,
         q50_c_cms = med_c * RI_10yr_Q_cms,
         q2.5_c_cms = q2.5_c * RI_10yr_Q_cms,
         q97.5_c_cms = q97.5_c * RI_10yr_Q_cms)

# Convert to cfs.
c_Q10_6 <- c_Q10_6 %>%
  mutate(mean_c_cfs = mean_c_cms * 35.3147,
         q50_c_cfs = q50_c_cms * 35.3147,
         q2.5_c_cfs = q2.5_c_cms * 35.3147,
         q97.5_c_cfs = q97.5_c_cms * 35.3147)

# Trim and export data.
c_Q10_6_trim <- c_Q10_6 %>%
  dplyr::select(site_name, mean_c_cfs, q2.5_c_cfs, q50_c_cfs, q97.5_c_cfs)

# write_csv(c_Q10_6_trim, 
#           "data_working/critical_discharge_cfs_6_sites_091222.csv")

#### Data Prep ####

# Re-running these 6 sites and normalizing by MAX Q rather than 10yr flood Q.
# Tried average Q, but the model freaked out because values were not between 0 and 1.

# First, filter available data by the sites of interest.
dat_in_4 <- dat_in %>%
  filter(site_name %in% my_list)

# And calculate mean/max Q for time periods for which we have data.
q_stats <- dat_in_4 %>%
  group_by(site_name) %>%
  summarize(meanQ = mean(Q),
            maxQ = max(Q))

# Join this info to the main dataset.
dat_in_4 <- left_join(dat_in_4, q_stats)

# And create the new columns with Q normalized to the mean & max.
dat_in_4 <- dat_in_4 %>%
  mutate(Q_rel_mean = Q/meanQ,
         Q_rel_max = Q/maxQ)

# Make as list and export for use on Teton.
dat4_list <- split(dat_in_4, dat_in_4$.id)

#saveRDS(dat4_list, "data_working/list_4sites_meanmaxQnorm_allSL.rds")

# And do the same for 2 additional sites per Jud's request.
# 0165389205 - Accotink Creek near Ranger Road, VA
# 05515500 - Kankakee River at Davis, IN

my_list2 <- c("nwis_0165389205", "nwis_05515500")

# First, filter available data by the 4 sites of interest.
dat_in_2 <- dat_in %>%
  filter(site_name %in% my_list2)

# And calculate mean Q for time periods for which we have data.
q_stats2 <- dat_in_2 %>%
  group_by(site_name) %>%
  summarize(meanQ = mean(Q),
            maxQ = max(Q))

# Join this info to the main dataset.
dat_in_2 <- left_join(dat_in_2, q_stats2)

# And create the new column with Q normalized to the mean.
dat_in_2 <- dat_in_2 %>%
  mutate(Q_rel_mean = Q/meanQ,
         Q_rel_max = Q/maxQ)

# Make as list and export for use on Teton.
dat2_list <- split(dat_in_2, dat_in_2$.id)

#saveRDS(dat2_list, "data_working/list_2sites_meanmaxQnorm_allSL.rds")

# And process outputs.
# Calculate normalized c statistic values.
# First, bind together data_out_max (4 sites) and data_out_max2 (2 sites).
dat_out_max_6 <- rbind(dat_out_max, dat_out_max2)

c_summary_max <- dat_out_max_6 %>%
  group_by(site_name) %>%
  summarize(mean_c = mean(c),
            med_c = median(c),
            q2.5_c = quantile(c, c(.025)),
            q97.5_c = quantile(c, c(.975))) %>%
  ungroup()

# Identify ten year flood values to back calculate c in cfs.
Qmax_4 <- dat_in2 %>%
  filter(site_name %in% my_list) %>%
  distinct(site_name,maxQ)

Qmax_2 <- dat_in3 %>%
  filter(site_name %in% my_list2) %>%
  distinct(site_name,maxQ)

Qmax_6 <- rbind(Qmax_4, Qmax_2)

# Join these datasets.
c_Qmax_6 <- left_join(c_summary_max, Qmax_6) # joined by "site_name"

# Calculate "c" in cubic meters per second.
c_Qmax_6 <- c_Qmax_6 %>%
  mutate(mean_c_cms = mean_c * maxQ,
         q50_c_cms = med_c * maxQ,
         q2.5_c_cms = q2.5_c * maxQ,
         q97.5_c_cms = q97.5_c * maxQ)

# Convert to cfs.
c_Qmax_6 <- c_Qmax_6 %>%
  mutate(mean_c_cfs = mean_c_cms * 35.3147,
         q50_c_cfs = q50_c_cms * 35.3147,
         q2.5_c_cfs = q2.5_c_cms * 35.3147,
         q97.5_c_cfs = q97.5_c_cms * 35.3147)

# Trim and export data.
c_Qmax_6_trim <- c_Qmax_6 %>%
  dplyr::select(site_name, mean_c_cfs, q2.5_c_cfs, q50_c_cfs, q97.5_c_cfs)

# write_csv(c_Qmax_6_trim, 
#           "data_working/critical_discharge_cfs_qmax_6_sites_091222.csv")

# Plot for lab meeting next week.
c_Qmax_6_trim$Normalization <- "Qmax"
c_Q10_6_trim$Normalization <- "Q10"

c_all <- rbind(c_Qmax_6_trim, c_Q10_6_trim)

(fig_c <- ggplot(c_all, aes(x = mean_c_cfs, y = site_name, color = Normalization)) +
    geom_point(size = 5, alpha = 0.75, position = position_dodge(width = -0.5)) +
    geom_errorbarh(aes(xmin = q2.5_c_cfs, xmax = q97.5_c_cfs), 
                   height = 0.25,
                   position = position_dodge(width = -0.5)) +
    labs(y = "Site", x = "c (cfs)") +
    scale_color_manual(values = c("#7AC9B7","#3793EC")) +
    theme_bw())

# export exploratory figures
# ggsave(("figures/teton_summer22/c_q10_vs_qmaxnorm_fig.png"),
#        width = 16,
#        height = 10,
#        units = "cm"
# )

# End of script.
