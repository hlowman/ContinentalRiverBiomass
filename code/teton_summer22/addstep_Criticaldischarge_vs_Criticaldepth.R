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
         "ggplot"), require, character.only=T)

#### Data Import ####

# Load site-level info (hypoxia and Appling datasets).
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
site <- fread("data_raw/site_data.tsv")

# Load site-level parameters dataset with all iterations included.
dat_out_q10 <- readRDS("data_working/teton_190rivers_model_site_params_all_iterations_082422.rds")
dat_out_max <- readRDS("data_working/teton_4rivers_model_site_params_all_iterations_090822.rds")

# Load in original data fed into the model.
dat_in <- readRDS("data_working/df_190sites_10yrQnorm_allSL.rds")
dat_in2 <- readRDS("data_working/list_4sites_meanmaxQnorm_allSL.rds")
dat_in2 <- map_df(dat_in2, ~as.data.frame(.x), .id="site_name")

#### Data Selection ####

# Geomorphic data is only available for the following sites:
# nwis_05451210: South Fork Iowa River, New Providence, IA
# nwis_05524500: Iroquois River, Foresman, IN
# nwis_05579630: Kickapoo Creek, Bloomington, IL
# nwis_01645704: Difficult Run, Fairfax, VA

# Filtering by those sites
my_list <- c("nwis_05451210","nwis_05524500",
             "nwis_05579630","nwis_01645704")

dat_out_q10_4 <- dat_out_q10 %>%
  filter(site_name %in% my_list)

# And calculate normalized c statistic values
c_summary_q10 <- dat_out_q10_4 %>%
  group_by(site_name) %>%
  summarize(mean_c = mean(c),
            med_c = median(c),
            q2.5_c = quantile(c, c(.025)),
            q97.5_c = quantile(c, c(.975))) %>%
  ungroup()

# Identify ten year flood values to back calculate c in cfs.
Q10_4 <- dat_in %>%
  filter(site_name %in% my_list) %>%
  distinct(site_name,RI_10yr_Q,RI_10yr_Q_cms)

# Join these datasets.
c_Q10_4 <- left_join(c_summary_q10, Q10_4) # joined by "site_name"

# Calculate "c" in cubic meters per second.
c_Q10_4 <- c_Q10_4 %>%
  mutate(mean_c_cms = mean_c * RI_10yr_Q_cms,
         q50_c_cms = med_c * RI_10yr_Q_cms,
         q2.5_c_cms = q2.5_c * RI_10yr_Q_cms,
         q97.5_c_cms = q97.5_c * RI_10yr_Q_cms)

# Convert to cfs.
c_Q10_4 <- c_Q10_4 %>%
  mutate(mean_c_cfs = mean_c_cms * 35.3147,
         q50_c_cfs = q50_c_cms * 35.3147,
         q2.5_c_cfs = q2.5_c_cms * 35.3147,
         q97.5_c_cfs = q97.5_c_cms * 35.3147)

# Trim and export data.
c_Q10_4_trim <- c_Q10_4 %>%
  dplyr::select(site_name, mean_c_cfs, q2.5_c_cfs, q50_c_cfs, q97.5_c_cfs)

write_csv(c_Q10_4_trim, 
          "data_working/critical_discharge_cfs_4_sites_090822.csv")

#### Data Prep ####

# Re-running these 4 sites and normalizing by MAX Q rather than 10yr flood Q.
# Tried average Q, but the model freaked out because values were not between 0 and 1.

# First, filter available data by the 4 sites of interest.
dat_in_4 <- dat_in %>%
  filter(site_name %in% my_list)

# And calculate mean Q for time periods for which we have data.
q_stats <- dat_in_4 %>%
  group_by(site_name) %>%
  summarize(meanQ = mean(Q),
            maxQ = max(Q))

# Join this info to the main dataset.
dat_in_4 <- left_join(dat_in_4, q_stats)

# And create the new column with Q normalized to the mean.
dat_in_4 <- dat_in_4 %>%
  mutate(Q_rel_mean = Q/meanQ,
         Q_rel_max = Q/maxQ)

# Make as list and export for use on Teton.
dat4_list <- split(dat_in_4, dat_in_4$.id)

#saveRDS(dat4_list, "data_working/list_4sites_meanmaxQnorm_allSL.rds")

# And process outputs.
# Calculate normalized c statistic values.

c_summary_max <- dat_out_max %>%
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

# Join these datasets.
c_Qmax_4 <- left_join(c_summary_max, Qmax_4) # joined by "site_name"

# Calculate "c" in cubic meters per second.
c_Qmax_4 <- c_Qmax_4 %>%
  mutate(mean_c_cms = mean_c * maxQ,
         q50_c_cms = med_c * maxQ,
         q2.5_c_cms = q2.5_c * maxQ,
         q97.5_c_cms = q97.5_c * maxQ)

# Convert to cfs.
c_Qmax_4 <- c_Qmax_4 %>%
  mutate(mean_c_cfs = mean_c_cms * 35.3147,
         q50_c_cfs = q50_c_cms * 35.3147,
         q2.5_c_cfs = q2.5_c_cms * 35.3147,
         q97.5_c_cfs = q97.5_c_cms * 35.3147)

# Trim and export data.
c_Qmax_4_trim <- c_Qmax_4 %>%
  dplyr::select(site_name, mean_c_cfs, q2.5_c_cfs, q50_c_cfs, q97.5_c_cfs)

write_csv(c_Qmax_4_trim, 
          "data_working/critical_discharge_cfs_qmax_4_sites_090822.csv")

# Plot for lab meeting next week.
c_Qmax_4_trim$Normalization <- "Qmax"
c_Q10_4_trim$Normalization <- "Q10"

c_all <- rbind(c_Qmax_4_trim, c_Q10_4_trim)

(fig_c <- ggplot(c_all, aes(x = mean_c_cfs, y = site_name, color = Normalization)) +
    geom_point(size = 5, alpha = 0.75, position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = q2.5_c_cfs, xmax = q97.5_c_cfs), 
                   height = 0.25,
                   position = position_dodge(width = 0.5)) +
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
