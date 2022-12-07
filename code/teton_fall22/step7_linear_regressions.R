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

# Next, the data for the Qc:Q2yr lm.
dat_Qc <- readRDS("data_working/QcQ2_filtered_141sites_113022.rds")

# Finally, the data for sites' Qc exceedances.
dat_exc <- readRDS("data_working/Qc_exceedances_141sites_120122.rds")

# And the hypoxia dataset for additional info re: precip. "pre_mm_cyr"
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")

# And the dataset with JUC10 delineations.
site_HUC <- readRDS("data_working/HUC12_159sites_120722.rds")

# Select for additional variables of interest.

# Annual average precip for local catchment/watershed from HydroATLAS (mm)
site_precip <- site_info %>%
  dplyr::select(SiteID, pre_mm_cyr, pre_mm_uyr)

site_HUC10 <- site_HUC %>%
  dplyr::select(site_name, huc10_id)

#### Model 1: rmax ####

# Trim imported data down to variables of interest.

dat_rmax_trim <- dat_rmax %>%
  dplyr::select(site_name,
         r_med:NHD_AREASQKM, NHD_RdDensCat:NHD_PctImp2011Ws, 
         Canal:Orthophosphate)

# Join with precip data.
dat_rmax_trim <- left_join(dat_rmax_trim, site_precip, by = c("site_name" = "SiteID"))

# Join with HUC10 data.
dat_rmax_trim <- left_join(dat_rmax_trim, site_HUC10)

# And visualize the relationships with rmax median values.

rmax_covs <- ggpairs(dat_rmax_trim %>% 
                       dplyr::select(-site_name,
                                     -huc10_id,
                                     -geometry))

# ggsave(rmax_covs,
#        filename = "figures/teton_fall22/rmax_covariates.jpg",
#        width = 50,
#        height = 50,
#        units = "cm")

# Some notes regarding these covariates.

# (1) I think I should log transform the following variables based on their
# relationships with r_med: NHD_AREASQKM, width_med, Nitrate, and orthoP.

# (2) Road density and impervious cover appear tightly correlated (0.838), so
# I should probably only include one of these in the final model build.

# (3) Then, summer temp and latitude appear tightly correlated (0.666).
# They do provide somewhat different information, but after a discussion
# with Joanna, I'll be removing lat and lon.

# (4) After that, NO3 and PO4 appear the most tightly correlated (0.548).
# Those aren't duplicates, so I will include them as an interaction term.

# (5) Remaining Pearson's correlation values are below 0.5.

# (6) Removing GPP because it's a value used to generate rmax, so therefore
# doesn't make sense to include in post-hoc analysis.

# Log transform necessary covariates.

dat_rmax_trim <- dat_rmax_trim %>%
  mutate(r_log = log10(r_med),
         #GPP_log = log10(meanGPP),
         area_log = log10(NHD_AREASQKM),
         width_log = log10(width_med),
         no3_log = log10(Nitrate),
         po4_log = log10(Orthophosphate)) %>%
  mutate(dam_num = as.numeric(as.character(Dam)))

# Notes on model structure:

# rmax ~ cvQ + light + temp + precip + size + roads + dams + 1 | HUC10

# I am going to include the following covariates as representatives of the
# corresponding environmental factors:

# cvQ - flow
# summer light - light at the stream surface during greatest canopy*
# daily light - light at the stream surface throughout the year*
# summer temp - water temperature (decoupled from air temperature)
# precipitation - aridity and land-water connectivity metric
# order - stream size**
# watershed area - stream size**
# width - stream size**
# road density - terrestrial development
# dam - aquatic development

# * Models will explore with summer vs. avg. daily light.
# ** Models will explore each of the size indicators separately.

# The following covariates have been removed for the following reasons:

# latitude/longitude - following discussion with Joanna, not a good measure
# GPP - a function of f(light, flow) and also a direct predictor of rmax
# % impervious land cover - too closely correlated with road density
# canal - another metric of terr/aq development but felt duplicative
# NO3/PO4 - 50% of sites have no data, so these will be examined separately

# One on one plots for covariates of interest vs. rmax.

plot(r_med ~ cvQ, data = dat_rmax_trim)
plot(r_med ~ meanL, data = dat_rmax_trim)
plot(r_med ~ summerL, data = dat_rmax_trim) # looks more related
plot(r_med ~ summerT, data = dat_rmax_trim)
plot(r_med ~ area_log, data = dat_rmax_trim)
plot(r_med ~ NHD_RdDensCat, data = dat_rmax_trim) # this too
plot(r_med ~ Dam, data = dat_rmax_trim)
plot(r_med ~ pre_mm_cyr, data = dat_rmax_trim)
hist(dat_rmax_trim$r_med)

# Going to log transform r_med too.
dat_rmax_trim <- dat_rmax_trim %>%
  mutate(logrmax = log10(r_med))

hist(dat_rmax_trim$logrmax)

# Ok, and making the final dataset with which to build models since I can't have
# NAs with many of these functions.
dat_rmax_lm <- dat_rmax_trim %>%
  dplyr::select(logrmax, cvQ, summerL, meanL, summerT, NHD_RdDensCat, Dam,
         Order, area_log, width_log, huc10_id) %>%
  drop_na() # 151 sites left

##### Step 1: Create lm() and check residuals.

# rmax ~ cvQ + light + temp + precip + size + roads + dams + 1 | HUC10

a1 <- lm(logrmax ~ cvQ + summerL + summerT + NHD_RdDensCat + Dam, 
         data = dat_rmax_lm)

plot(a1) # residuals looking better with the log transform

##### Step 2: Fit the lm() with GLS and compare to lme().

a2 <- gls(logrmax ~ cvQ + summerL + summerT + NHD_RdDensCat + Dam, 
          data = dat_rmax_lm) # effectively a lm

a3 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensCat + Dam, 
          random = ~1 | huc10_id, data = dat_rmax_lm) # with random effect

anova(a2, a3) # a3 preferred, and I want to keep the random term in

##### Step 3: Decide on a variance structure.

# Not going to add alternate variance structure at this time.

##### Step 4,5,6: Fit the lme(), compare with lm(), and check residuals.

# See Steps 2 & 3.

##### Step 7/8: Step-wise Optimal Fixed Structure.

a3.2 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensCat + Dam, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm) 

# Try alternate data for light.
a4 <- lme(logrmax ~ cvQ + meanL + summerT + NHD_RdDensCat + Dam, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm) 

anova(a3.2, a4) # a3.2 preferred - use summerL

# Try different covariates for stream size.
a5 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensCat + Dam + Order, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm) 

a6 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensCat + Dam + area_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm)

a7 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensCat + Dam + width_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm)

anova(a5, a6, a7) # order (a5) AIC is least, but see rationale below

# Moving forward with watershed area as the indicator of stream size.
# Since they're fairly similar results-wise, choose area, because it's
# the most objectively measured. Order can be a bit subjective/misleading
# and stream width was calculated, not measured directly..

# Investigate precipitation, but keep in mind, it drops 40+ records.
a6.2 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensCat + Dam + area_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_trim %>%
            dplyr::select(logrmax, cvQ, summerL, meanL, summerT, NHD_RdDensCat, 
                          Dam, Order, area_log, width_log, 
                          pre_mm_cyr, huc10_id) %>%
            drop_na())

a8 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensCat + Dam + area_log +
            pre_mm_cyr, 
            random = ~1 | huc10_id, 
            method = "ML",
            data = dat_rmax_trim %>%
            dplyr::select(logrmax, cvQ, summerL, meanL, summerT, NHD_RdDensCat, 
                          Dam, Order, area_log, width_log, 
                          pre_mm_cyr, huc10_id) %>%
              drop_na())

anova(a6.2, a8) # They are identical - precip adds nothing to the model.

# Investigate Qc exceedance metrics.

# Build a model to investigate Qc exceedance instead of cvQ as a measured of
# flow disturbance.

dat_rmax_trim2 <- left_join(dat_rmax_trim, dat_exc)

dat_rmax_lm2 <- dat_rmax_trim2 %>%
  dplyr::select(logrmax, cvQ, total_exc_events, total_exc_days, summerL, meanL, summerT, 
                NHD_RdDensCat, Dam, Order, area_log, width_log, huc10_id) %>%
  drop_na() # 134 sites left

a6.3 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensCat + Dam + area_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm2)

a9 <- lme(logrmax ~ total_exc_events + summerL + summerT + NHD_RdDensCat + 
            Dam + area_log, 
            random = ~1 | huc10_id, 
            method = "ML",
            data = dat_rmax_lm2)

a10 <- lme(logrmax ~ total_exc_days + summerL + summerT + NHD_RdDensCat + 
            Dam + area_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm2)

anova(a6.3, a9, a10) # exceedance events seems to be the best measure of discharge

summary(a9)

ggpairs(dat_rmax_lm2 %>%
          dplyr::select(logrmax, total_exc_events, summerL, summerT,
                 NHD_RdDensCat, Dam, area_log))

# Examine how the residuals model looks for the structure I'm working with.
# Fit linear model between rmax median values and cvQ.
fit1 <- lm(logrmax ~ cvQ, data = dat_rmax_lm)

# Add residuals to original dataset for fitting lm().
dat_rmax_lm$residuals <- residuals(fit1)

# Fit residuals and remove discharge metric.
a11 <- lme(residuals ~ summerL + summerT + NHD_RdDensCat + Dam + area_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm)

anova(a11) # No other metrics appear significant.

##### Step 9: Refit with REML for final full model

afinal <- lme(logrmax ~ total_exc_events + summerL + summerT + 
                NHD_RdDensCat + Dam + area_log, 
              random = ~1 | huc10_id, 
              method = "REML",
              data = dat_rmax_lm2)

# Examine model diagnostics.
plot(afinal, col = 1)
qqnorm(afinal)

# Examine contrasts between singular effects.
summary(afinal)

# Output model results.
anova(afinal)

#                  numDF denDF   F-value p-value
# (Intercept)          1   109 1056.4571  <.0001
# total_exc_events     1    16   19.5004  0.0004
# summerL              1    16    0.1382  0.7150
# summerT              1    16    0.4912  0.4935
# NHD_RdDensCat        1    16    0.9327  0.3485
# Dam                  3    16    0.1025  0.9574
# area_log             1    16    0.2410  0.6302

##### Nutrients #####

# And build separate model for nutrients.

dat_rmax_lm3 <- dat_rmax_trim2 %>%
  dplyr::select(logrmax, no3_log, po4_log, huc10_id) %>%
  drop_na() # 89 sites left

a12 <- lme(logrmax ~ no3_log + po4_log, 
           random = ~1 | huc10_id, 
           method = "ML",
           data = dat_rmax_lm3)

summary(a12)

# Refit with REML for final full model

afinal_nuts <- lme(logrmax ~ no3_log + po4_log,
              random = ~1 | huc10_id, 
              method = "REML",
              data = dat_rmax_lm3)

# Examine the output.
summary(afinal_nuts)

# Examine model diagnostics.
plot(afinal_nuts, col = 1)
qqnorm(afinal_nuts)

# Output model results.
anova(afinal_nuts)

#             numDF denDF  F-value p-value
# (Intercept)     1    73 863.2928  <.0001
# no3_log         1    13  12.6835  0.0035
# po4_log         1    13   0.4779  0.5016

#### Model 2: Qc:Q2yr ####

# Trim imported data down to variables of interest (17 total covariates).

dat_Qc_trim <- dat_Qc %>%
  select(site_name,
         c_med, Qc_Q2yr, cvQ, meanGPP:NHD_AREASQKM, 
         NHD_RdDensCat:NHD_PctImp2011Ws, 
         Canal:width_med)

# And visualize the relationships with Qc:Qy2r values.

QcQ2_covs <- ggpairs(dat_Qc_trim %>% select(-site_name) %>%
                       select(-NHD_RdDensWs, -NHD_PctImp2011Ws)) # trim for space

# ggsave(QcQ2_covs,
#        filename = "figures/teton_fall22/QcQ2_covariates.jpg",
#        width = 50,
#        height = 50,
#        units = "cm")

# Some notes regarding these covariates.

# (1) Road density and impervious cover appear tightly correlated (0.824), so
# I should probably only include one of these in the final model build again.

# (2) Stream width and meanGPP also appear correlated (0.618), more strongly
# than they did in the rmax covariate exploration.

# (3) Remaining Pearson's correlation values are below 0.5.

# Log transform necessary covariates.

dat_Qc_trim <- dat_Qc_trim %>%
  mutate(GPP_log = log10(meanGPP),
         area_log = log10(NHD_AREASQKM),
         width_log = log10(width_med))

# Notes on model structure:

# I am going to include the following covariates as representatives of the
# corresponding environmental factors:

# cvQ - flow
# GPP - biological productivity (function of flow & light)*
# latitude - a rough location/temperature index?
# longitude - a rough location/aridity index?
# order - stream size**
# watershed area - stream size**
# width - stream size**
# road density - terrestrial development
# dam - aquatic development

# * Models will explore with and without GPP since GPP ~ f(flow, light).
# ** Models will explore each of the size indicators separately.

# The following covariates have been removed for the following reasons:

# light - not a factor for flow disturbance thresholds
# temperature - not a factor for flow disturbance thresholds
# % impervious land cover - too closely correlated with road density
# canal - another metric of terr/aq development but felt duplicative
# NO3/PO4 - not a factor for flow disturbance thresholds

##### Stream Size #####

# Build initial set of models to investigate size covariates.

str(dat_Qc_trim) # keeping order and dam data as categorical

lm1_Qc <- lm(Qc_Q2yr ~ cvQ + GPP_log + Lat_WGS84 + Lon_WGS84 +
               Order + NHD_RdDensCat + Dam, 
               data = dat_Qc_trim)

lm2_Qc <- lm(Qc_Q2yr ~ cvQ + GPP_log + Lat_WGS84 + Lon_WGS84 +
               area_log + NHD_RdDensCat + Dam,  
               data = dat_Qc_trim)

lm3_Qc <- lm(Qc_Q2yr ~ cvQ + GPP_log + Lat_WGS84 + Lon_WGS84 +
               width_log + NHD_RdDensCat + Dam, 
               data = dat_Qc_trim)

# Examine the outputs.

summary(lm1_Qc)
summary(lm2_Qc)
summary(lm3_Qc)

# Examine the coefficients.

lm1_Qc_tidy <- broom::tidy(lm1_Qc) # using order
View(lm1_Qc_tidy) # Lon (p<0.05)

lm2_Qc_tidy <- broom::tidy(lm2_Qc) # using area
View(lm2_Qc_tidy) # Lon (p<0.05)

lm3_Qc_tidy <- broom::tidy(lm3_Qc) # using width
View(lm3_Qc_tidy) # Lon (p<0.05)

# Examine model fit.

lm1_Qc_fit <- broom::glance(lm1_Qc) # order
View(lm1_Qc_fit) # adj R2 = 0.03, sigma = 1.22, p = 0.23, nobs = 134

lm2_Qc_fit <- broom::glance(lm2_Qc) # area
View(lm2_Qc_fit) # adj R2 = 0.03, sigma = 1.22, p = 0.20, nobs = 134

lm3_Qc_fit <- broom::glance(lm3_Qc) # width
View(lm3_Qc_fit) #adj R2 = 0.03, sigma = 1.22, p = 0.17, nobs = 133

# Examine model diagnostics.
plot(lm1_Qc) # qqplot doesn't look good
plot(lm2_Qc) # looks the same
plot(lm3_Qc) # same outliers - 102, 122, 125

# Compare models.
aic1q <- AIC(lm1_Qc) # 450
aic2q <- AIC(lm2_Qc) # 446
aic3q <- AIC(lm3_Qc) # 443

# Moving forward with width as the indicator of stream size again.

# After examining the outliers of this model and finding they
# are all in Texas/arid West, I am going to try adding light back
# in the model.

# Need to pull summerL from rmax dataset.
dat_light <- dat_rmax_trim %>%
  select(site_name, summerL)

dat_Qc_trim_light <- left_join(dat_Qc_trim, dat_light)

# Build new model including light.

lm4_Qc <- lm(Qc_Q2yr ~ cvQ + GPP_log + Lat_WGS84 + Lon_WGS84 +
             summerL + width_log + NHD_RdDensCat + Dam, 
             data = dat_Qc_trim_light)

# Examine the outputs.

summary(lm4_Qc)

# Examine the coefficients.

lm4_Qc_tidy <- broom::tidy(lm4_Qc) # using width
View(lm4_Qc_tidy) # Lon (p<0.05)

# Examine model fit.

lm4_Qc_fit <- broom::glance(lm4_Qc) # order
View(lm4_Qc_fit) # adj R2 = 0.02, sigma = 1.23, p = 0.24, nobs = 133

# Examine model diagnostics.
plot(lm4_Qc) # diagnostic plots look the same, with the same outliers

# Compare models.
aic4q <- AIC(lm4_Qc) # 445 - so not really any better with light.

# So, moving forward with lm3_Qc as model of choice, although it did a 
# terrible job of explaining variance.

# Export table of current results.
write_csv(lm3_Qc_tidy, "data_working/lm_Qc.csv")

# End of script.