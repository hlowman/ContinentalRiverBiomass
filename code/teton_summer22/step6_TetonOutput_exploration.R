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
         "rstan","bayesplot","shinystan", "here",
         "GGally", "glmmTMB", "MuMIn", "effects",
         "DHARMa", "lme4", "multcomp", "patchwork",
         "calecopal"), require, character.only=T)

#### Data Import ####

# Load site-level info (hypoxia and Appling datasets).
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
site <- fread("data_raw/site_data.tsv")

# Load site-level parameters dataset with all iterations included.
dat <- readRDS("data_working/teton_190rivers_model_site_params_all_iterations_082422.rds")

# Load in original data fed into the model.
dat_in <- readRDS("data_working/df_190sites_10yrQnorm_allSL.rds")

#### Data Processing ####

# Create new dataset to pull out only r values.
dat_mean_r <- dat %>%
  group_by(site_name) %>%
  summarize(r_mean = mean(r, na.rm = TRUE)) %>%
  ungroup()

# Export this data for future use.
saveRDS(dat_mean_r, "data_working/teton_190rivers_mean_r_090122.rds")

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

# Also need to compile stats from the original dataset.
dat_in_summ <- dat_in %>%
  group_by(site_name) %>%
  summarize(PAR_surf_mean = mean(PAR_surface, na.rm = TRUE),
            Q_mean = mean(Q, na.rm = TRUE),
            CV_Q = (sd(Q, na.rm = TRUE)/mean(Q, na.rm = TRUE)))

# And join this with the larger dataset.
dat_exp <- left_join(dat_exp, dat_in_summ)

#### Exploration ####

# exploratory plots of rmax vs. lat, lon, elevation, dam, canal influences

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
  geom_point(color = "#97C2E2", alpha = 0.75) +
  labs(x = "Distance to Nearest Ditch/Canal\nas Quantile of Daily Average Distance\nto 80% Oxygen Turnover\n0 = Near, 95 = Far",
       y = "Maximum Growth Rate (rmax)") + 
  theme_bw())
# much lower influence of canals, but higher rmax with less canal influence

(fig6 <- ggplot(dat_exp, aes(x = struct.npdes_flag, y = r_mean)) +
    geom_point(color = "#69B9FA", alpha = 0.75) +
    labs(x = "Distance to Nearest NPDES\nas Quantile of Daily Average Distance\nto 80% Oxygen Turnover\n0 = Near, 95 = Far",
         y = "Maximum Growth Rate (rmax)") + 
    theme_bw())
# much lower influence of canals, but higher rmax with less canal influence

(fig7 <- ggplot(dat_exp, aes(x = NHD_AREASQKM, y = r_mean)) +
    geom_point(color = "#59A3F8", alpha = 0.75) +
    labs(x = "Watershed Area (sq km)",
         y = "Maximum Growth Rate (rmax)") + 
    theme_bw())
# rmax appears to increase with declining watershed area, but that may
# be due to clustering of smaller watersheds

(fig8 <- ggplot(dat_exp, aes(x = NHD_RdDensCat, y = r_mean)) +
    geom_point(color = "#4B8FF7", alpha = 0.75) +
    labs(x = "Road Density in Catchment",
         y = "Maximum Growth Rate (rmax)") + 
    theme_bw())
# rmax appears to increase with lower density of roads
# less roads = less flashy flows = higher rmax??

(fig9 <- ggplot(dat_exp, aes(x = LU_category, y = r_mean)) +
    geom_boxplot(color = "#5A7ECB", alpha = 0.75) +
    labs(x = "Land Use",
         y = "Maximum Growth Rate (rmax)") + 
    theme_bw())
# higher rmax in agricultural/urban, lower in forested/grassland

(fig9.2 <- ggplot(dat_exp, aes(x = LU_category, y = r_mean,
                               color = LU_category)) +
    geom_boxplot() +
    geom_jitter(alpha = 0.75, width = 0.2) +
    scale_color_manual(values = cal_palette("sierra1")) +
    labs(x = "Land Use",
         y = expression(Maximum~Growth~Rate~(r[max]))) + 
    theme_bw() +
    theme(legend.position = "none"))
# revised boxplot for lab meeting

# export exploratory figures
ggsave(("figures/teton_summer22/rmax_landuse_fig.png"),
       width = 12,
       height = 9,
       units = "cm"
)

(fig10 <- ggplot(dat_exp, aes(x = PAR_surf_mean, y = r_mean)) +
    geom_point(color = "#6B6D9F", alpha = 0.75) +
    labs(x = "Mean Daily PAR at Stream Surface",
         y = "Maximum Growth Rate (rmax)") + 
    theme_bw())
# a cluster at the bottom, but also a hump shape centered about 500?

(fig10.2 <- ggplot(dat_exp, aes(x = PAR_surf_mean, y = r_mean)) +
    geom_point(color = "#7AC9B7", alpha = 0.75, size = 3) +
    labs(x = "Mean Daily PAR at Stream Surface",
         y = expression(Maximum~Growth~Rate~(r[max]))) + 
    theme_bw())
# revised figure for lab meeting

(fig11 <- ggplot(dat_exp, aes(x = log10(Q_mean), y = r_mean)) +
    geom_point(color = "#4C4976", alpha = 0.75) +
    labs(x = "Log of Mean Daily Discharge",
         y = "Maximum Growth Rate (rmax)") + 
    theme_bw())
# with increasing discharge, increasing rmax?

(fig11.2 <- ggplot(dat_exp, aes(x = CV_Q, y = r_mean)) +
    geom_point(color = "#3793EC", alpha = 0.75, size = 3) +
    labs(x = "Coefficient of Variation of Discharge",
         y = expression(Maximum~Growth~Rate~(r[max])))+ 
    theme_bw())
# revised figure for lab meeting

(bernhardt_lab_fig <- fig10.2 + fig11.2)

# export exploratory figures
ggsave(("figures/teton_summer22/rmax_light_cvq_fig.png"),
       width = 18,
       height = 9,
       units = "cm"
)

(fig12 <- ggplot(dat_exp, aes(x = pre_mm_cyr, y = r_mean)) +
    geom_point(color = "#151E2F", alpha = 0.75) +
    labs(x = "Mean Annual Precipitation in the Catchment",
         y = "Maximum Growth Rate (rmax)") + 
    theme_bw())

fig_compiled <- fig12 + fig1 + fig2 + fig3 + fig4 + fig5 + fig6 +
  fig7 + fig8 + fig9 + fig10 + fig11 + plot_layout(ncol = 4)

# export exploratory figures
ggsave(("figures/teton_summer22/rmax_exploration_figs.png"),
       width = 32,
       height = 28,
       units = "cm"
)

#### Modeling ####

# Select for only columns of interest
covars <- dat_exp %>%
  dplyr::select(site_name, Lat_WGS84, Lon_WGS84, ele_mt_cav, NHD_AREASQKM,
         dis_m3_pyr, pre_mm_cyr, LU_category, NHD_RdDensCat,
         PAR_surf_mean, Q_mean)

#### Correlation ####

# check correlation of variables independent of rmax
(corr_sites <- ggpairs(covars %>% 
                         dplyr::select(-c(site_name, Lat_WGS84, 
                                          Lon_WGS84, ele_mt_cav))))

# discharge is correlated with a few other things

#### Linear Model ####

# log-transform prior to running model
dat_exp <- dat_exp %>%
  mutate(log_r = log10(r_mean),
         log_Q = log10(Q_mean)) %>%
  # and make land use a factor
  mutate(LU_cat_f = factor(LU_category))

# full model
m1 <- glmmTMB(log_r ~ NHD_AREASQKM + LU_cat_f + # subsidy vars
                NHD_RdDensCat + PAR_surf_mean + log_Q, #disturbance vars
              data = dat_exp)

# bare bones model
m2 <- glmmTMB(log_r ~ PAR_surf_mean + log_Q, #disturbance vars
              data = dat_exp)

AICc(m1, m2) # m1 better

# Output of the model.
# Note, summary() function looks at contrasts between singular effects.
summary(m1)

# discharge, road density, land use, and watershed size all significant

# Checking residuals. - need to figure out how in new package
#plot(m1, col = 1)
#qqnorm(m1)

# Post-hoc:
mHSD <- glht(m1, linfct=mcp(LU_cat_f="Tukey")) # Run a Tukey's post hoc analysis on land use.
summary(mHSD) 
# forested and grassland both have significantly higher rmax values
# than agricultural

# End of script.
