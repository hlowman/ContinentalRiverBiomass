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
         "modelsummary", "here", "nlme"), require, character.only=T)

#### Data ####

# Import necessary datasets.

# First, the data for the rmax lm.
dat_rmax <- readRDS("data_working/rmax_filtered_159sites_113022.rds")

# Next, the data for the Qc:Q2yr lm.
dat_Qc <- readRDS("data_working/QcQ2_filtered_141sites_113022.rds")

# Also , the data for maximum algal yields.
dat_yield <- readRDS("data_working/maxalgalyield_159sites_120822.rds")

# Finally, the data for sites' Qc exceedances.
dat_exc <- readRDS("data_working/Qc_exceedances_159sites_120822.rds")

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
                                     -huc10_id))

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
# Those aren't duplicates, so I will include them in separate models.

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
# longitude - as a geographic metric
# order - stream size**
# watershed area - stream size**
# width - stream size**
# road density - terrestrial development
# dam - aquatic development

# * Models will explore with summer vs. avg. daily light.
# ** Models will explore each of the size indicators separately.

# The following covariates have been removed for the following reasons:

# latitude - too closely correlated with temperature
# longitude - too much of a east coast data bias for this to be accurate
# GPP - a function of f(light, flow) and also a direct predictor of rmax
# % impervious land cover - too closely correlated with road density
# canal - another metric of terr/aq development but felt duplicative
# NO3/PO4 - 50% of sites have no data, so these will be examined separately

# One on one plots for covariates of interest vs. rmax.

hist(dat_rmax_trim$r_med)

# Going to log transform r_med too.
dat_rmax_trim <- dat_rmax_trim %>%
  mutate(logrmax = log10(r_med))

plot(logrmax ~ cvQ, data = dat_rmax_trim)
plot(logrmax ~ meanL, data = dat_rmax_trim)
plot(logrmax ~ summerL, data = dat_rmax_trim) # looks more related
plot(logrmax ~ summerT, data = dat_rmax_trim)
plot(logrmax ~ area_log, data = dat_rmax_trim)
plot(logrmax ~ Order, data = dat_rmax_trim)
plot(logrmax ~ NHD_RdDensWs, data = dat_rmax_trim) # this too
plot(logrmax ~ Dam, data = dat_rmax_trim)
plot(logrmax ~ pre_mm_cyr, data = dat_rmax_trim)
plot(logrmax ~ Lon_WGS84, data = dat_rmax_trim)
hist(dat_rmax_trim$logrmax)

# Ok, and making the final dataset with which to build models since I can't have
# NAs with many of these functions.
dat_rmax_lm <- dat_rmax_trim %>%
  dplyr::select(logrmax, cvQ, summerL, meanL, summerT, NHD_RdDensWs, 
                Dam, Order, area_log, width_log, huc10_id) %>%
  drop_na() # 151 sites left

##### Step 1: Create lm() and check residuals.

# rmax ~ cvQ + light + temp + longitude + size + roads + dams + 1 | HUC10

a1 <- lm(logrmax ~ cvQ + summerL + summerT + NHD_RdDensWs + Dam, 
         data = dat_rmax_lm)

plot(a1) # residuals looking better with the log transform

summary(a1) # examine initial model output without the grouping by watershed
# hmmm, let's see how this does with the other covariates added in

##### Step 2: Fit the lm() with GLS and compare to lme().

a2 <- gls(logrmax ~ cvQ + summerL + summerT + NHD_RdDensWs + Dam, 
          data = dat_rmax_lm) # effectively a lm

a3 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensWs + Dam, 
          random = ~1 | huc10_id, data = dat_rmax_lm) # with random effect

anova(a2, a3) # a3 preferred, and I want to keep the random term in

##### Step 3: Decide on a variance structure.

# After examining the standardized vs. fitted residuals for the below model
# structure, I decided to test out some variance structures.

# afinal <- lme(logrmax ~ total_exc_days + summerL + summerT + 
#                 area_log + RdDensWs_log + Dam,
#               random = ~1 | huc10_id, method = "REML", data = dat_rmax_lm2 %>%
#                 mutate(RdDensWs_log = log10(NHD_RdDensWs)))

# Model doesn't converge with longitude or light as variance effect.
# Residuals don't look any better with order or dam as variance effects.

a_wvar <- lme(logrmax ~ total_exc_days + summerL + summerT + 
                area_log + RdDensWs_log + Dam,
              random = ~1 | huc10_id, 
              weights = varIdent(form = ~1 | Dam),
            data = dat_rmax_lm2 %>%
                mutate(RdDensWs_log = log10(NHD_RdDensWs)))

a_wovar <- lme(logrmax ~ total_exc_days + summerL + summerT + 
                area_log + RdDensWs_log + Dam,
              random = ~1 | huc10_id,
              data = dat_rmax_lm2 %>%
                mutate(RdDensWs_log = log10(NHD_RdDensWs)))

anova(a_wvar, a_wovar)
plot(a_wvar)

# Explored some single-covariate model fits in the "one by one" section below,
# and ultimately landed on the fact that residuals look worse with less covariates.
# And the primary issue here may be the huge number of groupings (a.k.a. watersheds).

##### Step 4,5,6: Fit the lme(), compare with lm(), and check residuals.

# See Steps 2 & 3.

##### Step 7/8: Step-wise Optimal Fixed Structure.

a3.2 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensWs + Dam, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm) 

# Try alternate data for light.
a4 <- lme(logrmax ~ cvQ + meanL + summerT +NHD_RdDensWs + Dam, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm) 

anova(a3.2, a4) # a3.2 preferred - use summerL

# Try different covariates for stream size.
a5 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensWs + Dam + Order, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm) 

a6 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensWs + Dam + area_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm)

a7 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensWs + Dam + width_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm)

anova(a5, a6, a7) # order (a5) AIC is least, but moving forward with width
# since it's better than area and order is somewhat subjective.

# Investigate precipitation, but keep in mind, it drops 40+ records.
a6.2 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensWs + Dam + width_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_trim %>%
            dplyr::select(logrmax, cvQ, summerL, meanL, summerT, NHD_RdDensWs, 
                          Dam, Order, area_log, width_log, Lon_WGS84,
                          pre_mm_cyr, huc10_id) %>%
            drop_na())

a8 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensWs + Dam + width_log +
            pre_mm_cyr, 
            random = ~1 | huc10_id, 
            method = "ML",
            data = dat_rmax_trim %>%
            dplyr::select(logrmax, cvQ, summerL, meanL, summerT, NHD_RdDensWs, 
                          Dam, Order, area_log, width_log, Lon_WGS84,
                          pre_mm_cyr, huc10_id) %>%
              drop_na())

anova(a6.2, a8) # They are nearly identical, and precip is *somewhat* correlated with
# Lat/Lon, so since it forces us to drop so many records (n=115 remaining),
# I'm skipping it.

# Investigate Qc exceedance metrics.

# Build a model to investigate Qc exceedance instead of cvQ as a measured of
# flow disturbance.

dat_rmax_trim2 <- left_join(dat_rmax_trim, dat_exc)

dat_rmax_lm2 <- dat_rmax_trim2 %>%
  dplyr::select(logrmax, cvQ, total_exc_events, total_exc_days, summerL, summerT, 
                NHD_RdDensWs, Order, Dam, width_log, Lon_WGS84, huc10_id) %>%
  drop_na() # 151 sites left

a6.3 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensWs + 
              Dam + width_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm2)

a9 <- lme(logrmax ~ total_exc_events + summerL + summerT + NHD_RdDensWs + 
            Dam + width_log, 
            random = ~1 | huc10_id, 
            method = "ML",
            data = dat_rmax_lm2)

a10 <- lme(logrmax ~ total_exc_days + summerL + summerT + NHD_RdDensWs + 
            Dam + width_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm2)

anova(a6.3, a9, a10) # exceedance days seems to be slightly better measure of discharge

summary(a10)
summary(a6.3)

ggpairs(dat_rmax_lm2 %>%
          dplyr::select(logrmax, total_exc_days, summerL, summerT, Lon_WGS84,
                 NHD_RdDensWs, Dam, area_log))

# Examine how the residuals model looks for the structure I'm working with.
# Fit linear model between rmax median values and cvQ.
fit1 <- lm(logrmax ~ cvQ, data = dat_rmax_lm)

# Add residuals to original dataset for fitting lm().
dat_rmax_lm$residuals <- residuals(fit1)

# Fit residuals and remove discharge metric.
a11 <- lme(residuals ~ summerL + summerT + NHD_RdDensWs + Dam + width_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm)

summary(a11) # Doesn't appear to change all that much.

##### Step 9: Refit with REML for final full model

# First, need to examine log transformations for a few more covariates,
# after examining subsequent ggpairs() plots

a12 <- lme(logrmax ~ total_exc_days + summerL_log + summerT + 
             NHD_RdDensWs + Dam + width_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm2 %>%
            mutate(summerL_log = log10(summerL))) # log scaling light

a13 <- lme(logrmax ~ total_exc_days + summerL + summerT +  
             RdDensWs_log + Dam + width_log, 
           random = ~1 | huc10_id, 
           method = "ML",
           data = dat_rmax_lm2 %>%
             mutate(RdDensWs_log = log10(NHD_RdDensWs))) # log scaling roads

anova(a10, a12, a13) # Makes little difference either way.

# Precip didn't seem to add anything to the model, and subtracted >40 records,
# so choosing not to use it at this time. Nutrients examined below.
afinal <- lme(logrmax ~ total_exc_days + # rather than cvQ bc it uses Qc threshold
                summerL + # better performing than mean L
                summerT + # more informative than Lat, which it was correlated with
                width_log + # not best performing, but more trustworthy than order
                NHD_RdDensWs + # used by JB as metric of development, more sig. than Cat
                Dam, # more pertinent to our interests than "Canal"
              random = ~1 | huc10_id, # because ~30% of sites share a watershed
              method = "REML",
              data = dat_rmax_lm2) # n = 152

# Examine model diagnostics.
plot(afinal)
qqnorm(afinal)

# Examine model outputs.
aout <- summary(afinal)
coef(aout)

#                    Value    Std.Error  DF    t-value      p-value
# (Intercept)    -1.397494e+00 2.150842e-01 116 -6.4974274 2.128231e-09
# total_exc_days  1.141956e-03 3.173875e-04  26  3.5979875 1.321342e-03
# summerL         8.553462e-08 2.104997e-07  26  0.4063408 6.878175e-01
# summerT        -6.268022e-03 8.053115e-03  26 -0.7783351 4.433975e-01
# width_log       2.746842e-01 6.578702e-02  26  4.1753558 2.955439e-04
# NHD_RdDensWs    3.248718e-02 1.108714e-02  26  2.9301676 6.968060e-03
# Dam50          -9.908723e-02 9.266447e-02  26 -1.0693120 2.947565e-01
# Dam80          -2.693219e-02 8.595985e-02  26 -0.3133113 7.565438e-01
# Dam95           8.146615e-02 6.039462e-02  26  1.3488974 1.890018e-01

##### Additional one by one models #####

dat_rmax_lm3 <- dat_rmax_lm2 %>%
  mutate(RdDensWs_log = log10(NHD_RdDensWs))

a1.1 <- lme(logrmax ~ total_exc_events, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_lm3)

plot(a1.1) # WOW, that looks terrible. Maybe the zero-inflation here is the issue.

a1.2 <- lme(logrmax ~ cvQ, 
            random = ~1 | huc10_id, 
            method = "ML",
            data = dat_rmax_lm3)

plot(a1.2) # Same here though.

a1.3 <- lme(logrmax ~ summerL, 
            random = ~1 | huc10_id, 
            method = "ML",
            data = dat_rmax_lm3)

plot(a1.3) # Realizing this might be an issue of the random rather than fixed effects.

a1.4 <- lm(logrmax ~ summerL, data = dat_rmax_lm3)
plot(a1.4)

##### Nutrients #####

# And build separate model for nutrients.

dat_rmax_lm4 <- dat_rmax_trim2 %>%
  dplyr::select(logrmax, no3_log, po4_log, huc10_id) %>%
  drop_na() # 89 sites left

a14 <- lme(logrmax ~ no3_log + po4_log, 
           random = ~1 | huc10_id, 
           method = "ML",
           data = dat_rmax_lm4)

plot(a14) # Looks alright actually.
summary(a14)

# Refit with REML for final full model

afinal_nuts <- lme(logrmax ~ no3_log + po4_log,
              random = ~1 | huc10_id, 
              method = "REML",
              data = dat_rmax_lm4)

# Examine model diagnostics.
plot(afinal_nuts)
qqnorm(afinal_nuts) # Looking just fine.

# Examine model outputs.
aout_nuts <- summary(afinal_nuts)
coef(aout_nuts)

#                   Value  Std.Error DF   t-value      p-value
# (Intercept) -0.87892263 0.11910962 73 -7.379107 2.052277e-10
# no3_log      0.17520313 0.08306747 13  2.109167 5.489118e-02
# po4_log      0.05972981 0.08640591 13  0.691270 5.015545e-01

#### Model 2: Max. Algal Yield ####

# First, need to bind MAY estimates with remaining data
just_may <- dat_yield %>%
  dplyr::select(site_name, may_med)

dat_rmax_trim3 <- left_join(dat_rmax_trim2, just_may)

# And visualize the relationships with MAY median values.
may_covs <- ggpairs(dat_rmax_trim3 %>% 
                       dplyr::select(may_med, cvQ:Orthophosphate))

# ggsave(may_covs,
#        filename = "figures/teton_fall22/may_covariates.jpg",
#        width = 50,
#        height = 50,
#        units = "cm")

# Some notes regarding these covariates.

# (1) I think I will again log transform the following variables: 
# NHD_AREASQKM, width_med, Nitrate, and orthoP.

# (2) Will keep in mind correlations btw covariates found earlier.

# (3) CVq and stream width jump out as potentially important.

# Proposed starting model structure:

# max algal yield ~ cvQ + light + temp + size + roads + dams + 1 | HUC10

# One on one plots for covariates of interest vs. may.

hist(dat_rmax_trim3$may_med)

# Going to log transform may_med too.
dat_rmax_trim3 <- dat_rmax_trim3 %>%
  mutate(logmay = log10(may_med))

plot(logmay ~ cvQ, data = dat_rmax_trim3)
plot(logmay ~ summerL, data = dat_rmax_trim3)
plot(logmay ~ summerT, data = dat_rmax_trim3)
plot(logmay ~ width_log, data = dat_rmax_trim3)
plot(logmay ~ NHD_RdDensWs, data = dat_rmax_trim3)
plot(logmay ~ Dam, data = dat_rmax_trim3)
hist(dat_rmax_trim3$logmay)

# Ok, and making the final dataset with which to build models since I can't have
# NAs with many of these functions.
dat_may_lm <- dat_rmax_trim3 %>%
  dplyr::select(logmay, cvQ, summerL, summerT, NHD_RdDensWs, 
                Dam, width_log, huc10_id) %>%
  drop_na() # 151 sites left

##### Step 1: Create lm() and check residuals.

# max algal yield ~ cvQ + light + temp + size + roads + dams + 1 | HUC10

b1 <- lm(logmay ~ cvQ + summerL + summerT + width_log + NHD_RdDensWs + Dam, 
         data = dat_may_lm)

plot(b1) # residuals looking better with the log transform

summary(b1) # examine initial model output without the grouping by watershed
# width, temperature, and maybe roads/dams jump out as important?

##### Step 2: Fit the lm() with GLS and compare to lme().

b2 <- gls(logmay ~ cvQ + summerL + summerT + width_log + NHD_RdDensWs + Dam, 
          data = dat_may_lm) # effectively a lm

b3 <- lme(logmay ~ cvQ + summerL + summerT + width_log + NHD_RdDensWs + Dam, 
          random = ~1 | huc10_id, data = dat_may_lm) # with random effect

anova(b2, b3) # b3 preferred, and I want to keep the random term in

##### Step 3: Decide on a variance structure.

# None at the moment.

##### Step 4,5,6: Fit the lme(), compare with lm(), and check residuals.

# See Steps 2 & 3.

##### Step 7/8: Step-wise Optimal Fixed Structure.

# Investigate precipitation, but keep in mind, it drops 40+ records.
b4 <- lme(logmay ~ cvQ + summerL + summerT + width_log + NHD_RdDensWs + Dam, 
            random = ~1 | huc10_id, 
            method = "ML",
            data = dat_rmax_trim3 %>%
              dplyr::select(logmay, cvQ, summerL, summerT, NHD_RdDensWs, 
                            Dam, width_log, pre_mm_cyr, huc10_id) %>%
              drop_na())

b5 <- lme(logmay ~ cvQ + summerL + summerT + width_log + NHD_RdDensWs + Dam +
            pre_mm_cyr, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_rmax_trim3 %>%
            dplyr::select(logmay, cvQ, summerL, summerT, NHD_RdDensWs, 
                          Dam, width_log, pre_mm_cyr, huc10_id) %>%
            drop_na())

anova(b4, b5) # They're nearly identical, and precip is *somewhat* correlated with
# Lat/Lon, so since it forces us to drop so many records (n=115 remaining),
# I'm skipping it again.

# Investigate Qc exceedance metrics.

# Build a model to investigate Qc exceedance instead of cvQ as a measured of
# flow disturbance.

dat_rmax_trim4 <- left_join(dat_rmax_trim3, dat_exc)

dat_may_lm2 <- dat_rmax_trim4 %>%
  dplyr::select(logmay, cvQ, total_exc_events, total_exc_days, summerL, summerT, 
                NHD_RdDensWs, Dam, width_log, huc10_id) %>%
  drop_na() # 152 sites left

b6 <- lme(logmay ~ cvQ + summerL + summerT +  NHD_RdDensWs + Dam + width_log, 
            random = ~1 | huc10_id, 
            method = "ML",
            data = dat_may_lm2)

b7 <- lme(logmay ~ total_exc_events + summerL + summerT + NHD_RdDensWs + 
            Dam + width_log, 
          random = ~1 | huc10_id, 
          method = "ML",
          data = dat_may_lm2)

b8 <- lme(logmay ~ total_exc_days + summerL + summerT + NHD_RdDensWs + 
             Dam + width_log, 
           random = ~1 | huc10_id, 
           method = "ML",
           data = dat_may_lm2)

anova(b6, b7, b8) # exceedance events seems to be slightly better measure of discharge, going to include exceedance days for continuity with rmax

##### Step 9: Refit with REML for final full model

# Precip didn't seem to add anything to the model, and subtracted >40 records,
# so choosing not to use it at this time. Nutrients examined below.
bfinal <- lme(logmay ~ total_exc_days + # rather than cvQ bc it uses Qc threshold
                summerL + # better performing than mean L
                summerT + # more informative than Lat, which it was correlated with
                width_log + # not best performing, but most trustworthy size metric
                NHD_RdDensWs + # used by JB as metric of development
                Dam, # more pertinent to our interests than "Canal"
              random = ~1 | huc10_id, # because ~30% of sites share a watershed
              method = "REML",
              data = dat_may_lm2) # n = 151

# Examine model diagnostics.
plot(bfinal) # similar shape as rmax residuals
qqnorm(bfinal) # looks fine

# Examine model outputs.
bout <- summary(bfinal)
coef(bout)

#                        Value    Std.Error  DF    t-value      p-value
# (Intercept)    -1.461022e+00 4.444644e-01 116 -3.2871513 1.339547e-03
# total_exc_days  2.050021e-03 6.559679e-04  26  3.1251851 4.333923e-03
# summerL         4.917369e-07 4.355078e-07  26  1.1291116 2.691612e-01
# summerT        -3.292882e-02 1.662697e-02  26 -1.9804455 5.832675e-02
# width_log       6.928018e-01 1.359675e-01  26  5.0953496 2.613403e-05
# NHD_RdDensWs    4.451154e-02 2.293182e-02  26  1.9410382 6.316716e-02
# Dam50          -7.789439e-02 1.928611e-01  26 -0.4038886 6.895975e-01
# Dam80           1.292574e-01 1.781592e-01  26  0.7255162 4.746156e-01
# Dam95           2.898039e-01 1.251119e-01  26  2.3163581 2.868362e-02

##### Nutrients #####

# And build separate model for nutrients.

dat_may_lm3 <- dat_rmax_trim3 %>%
  dplyr::select(logmay, no3_log, po4_log, huc10_id) %>%
  drop_na() # 89 sites left

b9 <- lme(logmay ~ no3_log + po4_log, 
           random = ~1 | huc10_id, 
           method = "ML",
           data = dat_may_lm3)

plot(b9) # Looks ok.
summary(b9)

# Refit with REML for final full model

bfinal_nuts <- lme(logmay ~ no3_log + po4_log,
                   random = ~1 | huc10_id, 
                   method = "REML",
                   data = dat_may_lm3)

# Examine model diagnostics.
plot(bfinal_nuts)
qqnorm(bfinal_nuts) # Looking just fine.

# Examine model outputs.
bout_nuts <- summary(bfinal_nuts)
coef(bout_nuts)

#                   Value Std.Error DF    t-value     p-value
# (Intercept) -0.73426230 0.2560306 73 -2.8678691 0.005398894
# no3_log      0.32742324 0.1789020 13  1.8301825 0.090238344
# po4_log      0.09956532 0.1856118 13  0.5364169 0.600729814


#### Model 3: Qc:Q2yr ####

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