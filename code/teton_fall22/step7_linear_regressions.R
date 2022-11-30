## Resilience of Stream Productivity to Disturbance
## November 30, 2022
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

# This file will perform the linear regressions to explore covariates that
# might explain the median parameter estimates from our model output.

#### Setup ####

# Load necessary packages.
lapply(c("tidyverse", "lubridate", "data.table",
         "GGally", "glmmTMB", "MuMIn", "effects",
         "DHARMa", "lme4", "multcomp", "patchwork",
         "modelsummary", "here"), require, character.only=T)

#### Data ####

# Import necessary datasets.

# First, the data for the rmax lm.
dat_rmax <- readRDS("data_working/rmax_filtered_159sites_113022.rds")

#### Model 1: rmax ####

# Trim imported data down to variables of interest (17 total covariates).

dat_rmax_trim <- dat_rmax %>%
  select(site_name,
         r_med:NHD_AREASQKM, NHD_RdDensCat:NHD_PctImp2011Ws, 
         Canal:Orthophosphate)

# And visualize the relationships with rmax median values.

rmax_covs <- ggpairs(dat_rmax_trim %>% select(-site_name) %>%
          select(-meanL, -NHD_RdDensWs, -NHD_PctImp2011Ws)) # trimming for space

# ggsave(rmax_covs,
#        filename = "figures/teton_fall22/rmax_covariates.jpg",
#        width = 50,
#        height = 50,
#        units = "cm")

# Some notes regarding these covariates.

# (1) I think I should log transform the following variables based on their
# relationships with r_med: meanGPP, NHD_AREASQKM, width_med, Nitrate,
# and Orthophosphate.

# (2) Road density and impervious cover appear tightly correlated (0.838), so
# I should probably only include one of these in the final model build.

# (3) Then, summer temperature and latitude appear tightly correlated (0.666).
# They do provide somewhat different information, so I will include them as
# an interaction term. (Thinking about latitude as a temp gradient while
# longitude might be a kind of aridity gradient.)

# (4) After that, nitrate and orthoP appear the most tightly correlated (0.548).
# Those aren't duplicates, so I will include them as an interaction term.

# (5) Remaining Pearson's correlation values are below 0.5.

# Log transform necessary covariates.

dat_rmax_trim <- dat_rmax_trim %>%
  mutate(GPP_log = log10(meanGPP),
         area_log = log10(NHD_AREASQKM),
         width_log = log10(width_med),
         no3_log = log10(Nitrate),
         po4_log = log10(Orthophosphate))

# Build initial, saturated model.
# The nutrient data was causing nearly half of the observations (72)
# to be dropped due to complete case analysis. So, these will
# be examined in another separate model.

str(dat_rmax_trim) # keeping order, canal, and dam data as categorical

lm1_rmax <- lm(r_med ~ cvQ + GPP_log + summerL + summerT*Lat_WGS84 +
                 Lon_WGS84 + Order + area_log + NHD_RdDensCat + 
                 Canal + Dam + width_log, data = dat_rmax_trim)

# Examine the outputs.

summary(lm1_rmax)

# Examine the coefficients.

lm1_rmax_tidy <- broom::tidy(lm1_rmax)
View(lm1_rmax_tidy) # GPP, order (4,5,6,9), longitude, and dam (80%)
# appear to be significant factors

# Examine model fit.

lm1_rmax_fit <- broom::glance(lm1_rmax)
View(lm1_rmax_fit) # nobs = 151

# Examine model diagnostics.
plot(lm1_rmax) # site 30 does appear to be an outlier here

# Build trimmed down model.
# Removing Canal because I don't feel it needs to be in here.
# Removing Lat because it was closely coupled with temperature.
# Removing area_log because it's close to stream size in concept.
# Also trying without order but keeping in width_log for stream size.

lm2_rmax <- lm(r_med ~ cvQ + GPP_log + summerL + summerT +
                 Lon_WGS84 + NHD_RdDensCat + Dam + Order, 
               data = dat_rmax_trim)

lm3_rmax <- lm(r_med ~ cvQ + GPP_log + summerL + summerT +
                 Lon_WGS84 + NHD_RdDensCat + Dam + width_log, 
               data = dat_rmax_trim)

# Examine the outputs.

summary(lm2_rmax)
summary(lm3_rmax)

# And build separate model for nutrients.
lmnuts_rmax <- lm(r_med ~ no3_log*po4_log, 
                  data = dat_rmax_trim)

# Examin the outputs.
summary(lmnuts_rmax)

# Examine the coefficients.

lmnuts_rmax_tidy <- broom::tidy(lmnuts_rmax)
View(lmnuts_rmax_tidy) # neither NO3 nor PO4 is significant

# Examine model fit.

lmnuts_rmax_fit <- broom::glance(lmnuts_rmax)
View(lmnuts_rmax_fit) # adj R2 = 0.09, sigma = 0.09, p = 0.01, nobs = 89

# Examine model diagnostics.
plot(lmnuts_rmax) # site 135 does *nearly* appear to be an outlier

# End of script.