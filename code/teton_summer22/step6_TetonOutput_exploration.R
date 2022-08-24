## Resilience of Stream Productivity to Disturbance
## August 22, 2022
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

#### Setup ####

## Load packages
lapply(c("tidyverse", "lubridate", "data.table",
         "rstan","bayesplot","shinystan", "here"), require, character.only=T)

#### Data Import ####

# Load site-level info (hypoxia and Appling datasets).
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
site <- fread("data_raw/site_data.tsv")

# Load site-level parameters dataset with all iterations included.
dat <- readRDS("data_working/teton_190rivers_model_site_params_all_iterations_082422.rds")

#### Data Processing ####

# Create new dataset to pull out only r values.
dat_mean_r <- dat %>%
  group_by(site_name) %>%
  summarize(r_mean = mean(r, na.rm = TRUE)) %>%
  ungroup()

# And remove 10 sites that didn't pass the diagnostics checks:
to_remove <- c("nwis_0165389205", "nwis_02171645", "nwis_02336300",
               "nwis_02336728", "nwis_03067510", "nwis_03073000",
               "nwis_03289200", "nwis_03302030", "nwis_04124000",
               "nwis_06893350" )

dat_mean_r180 <- dat_mean_r %>%
  filter(!site_name %in% to_remove) %>%
  # as well as any values below 0 (15 sites)
  filter(r_mean > 0)

# And join the datasets together.
dat_exp <- left_join(dat_mean_r180, site_info, by = c("site_name"="SiteID"))
dat_exp <- left_join(dat_exp, site)

#### Exploration ####

# exploratory plots of rmax vs. lat, lon, elevation, dam, and canal influences

(fig1 <- ggplot(dat_exp, aes(x = Lat_WGS84, y = r_mean)) +
  geom_point(color = "#E29244", alpha = 0.75) +
  labs(x = "Latitude",
       y = "Maximum Growth Rate (rmax)") + 
  theme_bw())
# rmax increases at higher latitudes (north)

(fig2 <- ggplot(dat_exp, aes(x = Lon_WGS84, y = r_mean)) +
  geom_point(color = "#FFAA00", alpha = 0.75) +
  labs(x = "Longitude",
       y = "Maximum Growth Rate (rmax)") + 
  theme_bw())
# rmax increases at lower longitudes (east)

(fig3 <- ggplot(dat_exp, aes(x = ele_mt_cav, y = r_mean)) +
  geom_point(color = "#D46F10", alpha = 0.75) +
  labs(x = "Mean Catchment Elevation (m)",
       y = "Maximum Growth Rate (rmax)") + 
  theme_bw())
# clusters at low elevations, so distinct pattern not clear
# 40 sites missing elevation data

(fig4 <- ggplot(dat_exp, aes(x = struct.dam_flag, y = r_mean)) +
  geom_point(color = "#4CA49E", alpha = 0.75) +
  labs(x = "Distance to Nearest Dam\nas Quantile of Daily Average Distance\nto 80% Oxygen Turnover\n0 = Near, 95 = Far",
       y = "Maximum Growth Rate (rmax)") + 
  theme_bw())
# unclear if there is any pattern w/ dams (which is still interesting!)

(fig5 <- ggplot(dat_exp, aes(x = struct.canal_flag, y = r_mean)) +
  geom_point(color = "#69B9FA", alpha = 0.75) +
  labs(x = "Distance to Nearest Ditch/Canal\nas Quantile of Daily Average Distance\nto 80% Oxygen Turnover\n0 = Near, 95 = Far",
       y = "Maximum Growth Rate (rmax)") + 
  theme_bw())
# much lower influence of canals, but higher rmax with less canal influence

(fig6 <- ggplot(dat_exp, aes(x = NHD_AREASQKM, y = r_mean)) +
    geom_point(color = "#59A3F8", alpha = 0.75) +
    labs(x = "Watershed Area (sq km)",
         y = "Maximum Growth Rate (rmax)") + 
    theme_bw())
# rmax appears to increase with declining watershed area, but that may
# be due to clustering of smaller watersheds

(fig7 <- ggplot(dat_exp, aes(x = NHD_RdDensCat, y = r_mean)) +
    geom_point(color = "#4B8FF7", alpha = 0.75) +
    labs(x = "Road Density in Catchment",
         y = "Maximum Growth Rate (rmax)") + 
    theme_bw())
# rmax appears to increase with lower density of roads
# less roads = less flashy flows = higher rmax??

(fig8 <- ggplot(dat_exp, aes(x = LU_category, y = r_mean)) +
    geom_boxplot(color = "#5A7ECB", alpha = 0.75) +
    labs(x = "Land Use",
         y = "Maximum Growth Rate (rmax)") + 
    theme_bw())
# higher rmax in agricultural/urban, lower in forested/grassland

fig_compiled <- (fig1 + fig2) / (fig3 + fig4) / (fig5 + fig6) / (fig7 + fig8)

# export exploratory figures
# ggsave(("figures/teton_summer22/rmax_exploration_figs.png"),
#        width = 20,
#        height = 32,
#        units = "cm"
# )

# End of script.
