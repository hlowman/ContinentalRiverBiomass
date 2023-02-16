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

# This file will perform the post-hoc analyses to explore covariates that
# might explain the median parameter estimates from our model output.

#### Setup ####

# Load necessary packages.
lapply(c("tidyverse", "lubridate", "data.table",
         "GGally", "glmmTMB", "MuMIn", "effects",
         "DHARMa", "lme4", "multcomp", "patchwork",
         "modelsummary", "here", "nlme", "brms", "loo",
         "bayesplot", "tidybayes"), require, character.only=T)

#### Data ####

# Import necessary datasets.

# First, the data for the rmax lm.
dat_rmax <- readRDS("data_working/rmax_filtered_159sites_113022.rds")

# Next, the data for the Qc:Q2yr lm.
dat_Qc <- readRDS("data_working/QcQ2_filtered_141sites_113022.rds")

# Also , the data for maximum algal yields.
dat_yield <- readRDS("data_working/maxalgalyield_159sites_021323.rds")

# Finally, the data for sites' Qc exceedances.
dat_exc <- readRDS("data_working/Qc_exceedances_159sites_021423.rds")

# And the hypoxia dataset for additional info re: precip. "pre_mm_cyr"
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")

# And the dataset with all HUC delineations.
site_HUC <- readRDS("data_working/HUC12_159sites_120922.rds")

# Select for additional variables of interest.

# Annual average precip for local catchment/watershed from HydroATLAS (mm)
site_precip <- site_info %>%
  dplyr::select(SiteID, pre_mm_cyr, pre_mm_uyr)

site_HUC2 <- site_HUC %>%
  dplyr::select(site_name, huc2_id)

#### Model 1: rmax LMEM  ####

# Trim imported data down to variables of interest.

dat_rmax_trim <- dat_rmax %>%
  dplyr::select(site_name,
         r_med:NHD_AREASQKM, NHD_RdDensCat:NHD_PctImp2011Ws, 
         Canal:Orthophosphate)

# Join with precip data.
dat_rmax_trim <- left_join(dat_rmax_trim, site_precip, by = c("site_name" = "SiteID"))

# Join with HUC2 data.
dat_rmax_trim <- left_join(dat_rmax_trim, site_HUC2)

# And visualize the relationships with rmax median values.

rmax_covs <- ggpairs(dat_rmax_trim %>% 
                       dplyr::select(-site_name,
                                     -huc2_id))

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

# (6) Removing GPP, light, and Q, because they were used to generate rmax,
# so therefore doesn't make sense to include in post-hoc analysis.

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

# rmax ~ temp + precip + size + roads + dams + nuts + 1 | HUC2

# I am going to include the following covariates as representatives of the
# corresponding environmental factors:

# summer temp - water temperature (decoupled from air temperature)
# precipitation - aridity and land-water connectivity metric
# order - stream size**
# watershed area - stream size**
# width - stream size**
# road density - terrestrial development
# dam - aquatic development
# nutrients - NO3 + PO4 where available
# HUC2 - to account for random effect of geography

# ** Models will explore each of the size indicators separately.

# The following covariates have been removed for the following reasons:

# latitude - too closely correlated with temperature
# longitude - too much of a east coast data bias for this to be accurate
# GPP - a function of f(light, flow) and also a direct predictor of rmax
# cvQ - flow but also a direct predictor of rmax
# summer light - light at the stream surface during greatest canopy, but
# also a direct predictor of biomass/rmax
# daily light - light at the stream surface throughout the year, but
# also a direct predictor of biomass/rmax
# % impervious land cover - too closely correlated with road density
# canal - another metric of terr/aq development but felt duplicative

# One on one plots for covariates of interest vs. rmax.

hist(dat_rmax_trim$r_med)

# Going to log transform r_med too.
dat_rmax_trim <- dat_rmax_trim %>%
  mutate(logrmax = log10(r_med))

#plot(logrmax ~ cvQ, data = dat_rmax_trim)
#plot(logrmax ~ meanL, data = dat_rmax_trim)
#plot(logrmax ~ summerL, data = dat_rmax_trim)
plot(logrmax ~ summerT, data = dat_rmax_trim)
plot(logrmax ~ area_log, data = dat_rmax_trim)
plot(logrmax ~ Order, data = dat_rmax_trim)
plot(logrmax ~ width_log, data = dat_rmax_trim)
plot(logrmax ~ NHD_RdDensWs, data = dat_rmax_trim)
plot(logrmax ~ Dam, data = dat_rmax_trim)
plot(logrmax ~ pre_mm_cyr, data = dat_rmax_trim)
plot(logrmax ~ huc2_id, data = dat_rmax_trim)
hist(dat_rmax_trim$logrmax)

# Ok, and making the final dataset with which to build models since I can't have
# NAs with many of these functions.
dat_rmax_lm <- dat_rmax_trim %>%
  dplyr::select(logrmax, summerT, Order, area_log, width_log, NHD_RdDensWs, 
                Dam, huc2_id, cvQ) %>%
  drop_na() # 151 sites left

##### Step 1: Create lm() and check residuals.

# rmax ~ temp + roads + dams

a1 <- lm(logrmax ~ summerT + NHD_RdDensWs + Dam, 
         data = dat_rmax_lm)

plot(a1) # residuals looking better with the log transform

summary(a1) # examine initial model output without the grouping by watershed
# hmmm, let's see how this does with the other covariates added in

##### Step 2: Fit the lm() with GLS and compare to lme().

a2 <- gls(logrmax ~ summerT + NHD_RdDensWs + Dam, 
          data = dat_rmax_lm) # effectively a lm

a3 <- lme(logrmax ~ summerT + NHD_RdDensWs + Dam, 
          random = ~1 | huc2_id, data = dat_rmax_lm) # with random effect

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

# a_wvar <- lme(logrmax ~ total_exc_days + summerL + summerT + 
#                 area_log + RdDensWs_log + Dam,
#               random = ~1 | huc10_id, 
#               weights = varIdent(form = ~1 | Dam),
#             data = dat_rmax_lm2 %>%
#                 mutate(RdDensWs_log = log10(NHD_RdDensWs)))

# a_wovar <- lme(logrmax ~ total_exc_days + summerL + summerT + 
#                 area_log + RdDensWs_log + Dam,
#               random = ~1 | huc10_id,
#               data = dat_rmax_lm2 %>%
#                 mutate(RdDensWs_log = log10(NHD_RdDensWs)))

# anova(a_wvar, a_wovar)
# plot(a_wvar)

# Explored some single-covariate model fits in the "one by one" section below.
# Ultimately landed on the fact that residuals look worse with less covariates.
# And the primary issue here may be the huge number of groupings (a.k.a. 
# watersheds). So, went back and got larger HUCs to group by.

##### Step 4,5,6: Fit the lme(), compare with lm(), and check residuals.

# See Steps 2 & 3.

##### Step 7/8: Step-wise Optimal Fixed Structure.

# Try different covariates for stream size.
a4 <- lme(logrmax ~ summerT + NHD_RdDensWs + Dam + Order, 
          random = ~1 | huc2_id, 
          method = "ML",
          data = dat_rmax_lm) # stream order

a5 <- lme(logrmax ~ summerT + NHD_RdDensWs + Dam + area_log, 
          random = ~1 | huc2_id, 
          method = "ML",
          data = dat_rmax_lm) # watershed area

a6 <- lme(logrmax ~ summerT + NHD_RdDensWs + Dam + width_log, 
          random = ~1 | huc2_id, 
          method = "ML",
          data = dat_rmax_lm) # stream width

anova(a4, a5, a6) # order (a5) AIC is least, but moving forward with width
# since it's better than area and order is somewhat subjective.

# Investigate precipitation, but keep in mind, it drops 40+ records. (Also,
# chose to archive this because investigating precip separately below.)
# a6.2 <- lme(logrmax ~ summerT + NHD_RdDensWs + Dam + width_log, 
#           random = ~1 | huc2_id, 
#           method = "ML",
#           data = dat_rmax_trim %>%
#             dplyr::select(logrmax, summerT, NHD_RdDensWs, 
#                           Dam, width_log, pre_mm_cyr, huc2_id) %>%
#             drop_na())

# a7 <- lme(logrmax ~ summerT + NHD_RdDensWs + Dam + width_log +
#             pre_mm_cyr, 
#             random = ~1 | huc2_id, 
#             method = "ML",
#             data = dat_rmax_trim %>%
#             dplyr::select(logrmax, summerT, NHD_RdDensWs, 
#                           Dam, width_log, pre_mm_cyr, huc2_id) %>%
#               drop_na())

# anova(a6.2, a7) # Nearly identical, and precip is *somewhat* correlated with
# Lat/Lon, so since it forces us to drop so many records (n=115 remaining),
# I'm skipping it.

# Previous investigation Qc exceedance metrics. (Chose to archive this since
# I'm not putting a measure of discharge in the post-hoc model.)

# Build a model to investigate Qc exceedance instead of cvQ as a measured of
# flow disturbance.

#dat_rmax_trim2 <- left_join(dat_rmax_trim, dat_exc)

# dat_rmax_lm2 <- dat_rmax_trim2 %>%
#   dplyr::select(logrmax, cvQ, total_exc_events, total_exc_days, summerL, summerT, 
#                 NHD_RdDensWs, Order, Dam, width_log, Lon_WGS84, huc2_id) %>%
#   drop_na() # 151 sites left

# a6.3 <- lme(logrmax ~ cvQ + summerL + summerT + NHD_RdDensWs + 
#               Dam + width_log, 
#           random = ~1 | huc2_id, 
#           method = "ML",
#           data = dat_rmax_lm2)

# a9 <- lme(logrmax ~ total_exc_events + summerL + summerT + NHD_RdDensWs + 
#             Dam + width_log, 
#             random = ~1 | huc2_id, 
#             method = "ML",
#             data = dat_rmax_lm2)

# a10 <- lme(logrmax ~ total_exc_days + summerL + summerT + NHD_RdDensWs + 
#             Dam + width_log, 
#           random = ~1 | huc2_id, 
#           method = "ML",
#           data = dat_rmax_lm2)

#anova(a6.3, a9, a10) # exceedance events seems to be slightly better measure of discharge

#summary(a9)
#summary(a6.3)

# Examine how the residuals model looks for the structure I'm working with.
# Fit linear model between rmax median values and cvQ.
fit1 <- lm(logrmax ~ cvQ, data = dat_rmax_lm)

# Add residuals to original dataset for fitting lm().
dat_rmax_lm$residuals <- residuals(fit1)

# Fit residuals and remove discharge metric.
a8 <- lme(residuals ~ summerT + NHD_RdDensWs + Dam + width_log, 
          random = ~1 | huc2_id, 
          method = "ML",
          data = dat_rmax_lm)

summary(a6)
summary(a8) # Doesn't appear to change all that much.

##### Step 9: Refit with REML for final full model

# Maximize data availability:
dat_rmax_lm2 <- dat_rmax_trim %>%
  dplyr::select(logrmax, summerT, width_log, NHD_RdDensWs, Dam, huc2_id) %>%
  drop_na() #151

# Precip didn't seem to add anything to the model, and subtracted >40 records,
# so choosing not to use it at this time. Nutrients examined below.
afinal <- lme(logrmax ~ # log-transformed because right-skewed
                summerT + # more informative than Lat, which it was correlated with
                NHD_RdDensWs + # used by JB as metric of development, more sig. than Cat
                Dam + # more pertinent to our interests than "Canal"
                width_log, # not best performing, but more trustworthy than order
              random = ~1 | huc2_id, # because HUC10 was too few sites/watershed
              method = "REML",
              data = dat_rmax_lm) # n = 151

# Examine model diagnostics.
plot(afinal) # YAY! Residuals looking WAY better using HUC2.
qqnorm(afinal) # And qqplot looks just fine.

# Examine model outputs.
aout <- summary(afinal)
coef(aout)

#                    Value   Std.Error  DF    t-value      p-value
# (Intercept)  -1.17660341 0.237187782 130 -4.9606409 2.155395e-06
# summerT      -0.01774056 0.009215881 130 -1.9249985 5.641388e-02
# NHD_RdDensWs  0.04005887 0.011634177 130  3.4432064 7.739853e-04
# Dam50        -0.15159227 0.101653185 130 -1.4912693 1.383134e-01
# Dam80        -0.05272356 0.087398021 130 -0.6032581 5.473881e-01
# Dam95         0.04216065 0.062587345 130  0.6736289 5.017432e-01
# width_log     0.33187653 0.071093583 130  4.6681643 7.464196e-06

r.squaredGLMM(afinal)

#           R2m       R2c
# [1,] 0.171693 0.3000676

##### Additional one by one models #####

# These were done when I was trying to figure out what the issue was with
# poor residual distribution.

# dat_rmax_lm3 <- dat_rmax_lm2 %>%
#   mutate(RdDensWs_log = log10(NHD_RdDensWs))
# 
# a1.1 <- lme(logrmax ~ total_exc_events, 
#           random = ~1 | huc10_id, 
#           method = "ML",
#           data = dat_rmax_lm3)
# 
# plot(a1.1) # WOW, that looks terrible. Maybe the zero-inflation here is the issue.
# 
# a1.2 <- lme(logrmax ~ cvQ, 
#             random = ~1 | huc10_id, 
#             method = "ML",
#             data = dat_rmax_lm3)
# 
# plot(a1.2) # Same here though.
# 
# a1.3 <- lme(logrmax ~ summerL, 
#             random = ~1 | huc10_id, 
#             method = "ML",
#             data = dat_rmax_lm3)
# 
# plot(a1.3) # Realizing this might be an issue of the random rather than fixed effects.
# 
# a1.4 <- lm(logrmax ~ summerL, data = dat_rmax_lm3)
# plot(a1.4)

##### Nutrients #####

# And build separate model for nutrients.

dat_rmax_lm3 <- dat_rmax_trim %>%
  dplyr::select(logrmax, summerT, width_log, NHD_RdDensWs, 
                Dam, huc2_id, no3_log, po4_log) %>%
  drop_na() # 87 sites left

a9 <- lme(logrmax ~ summerT + NHD_RdDensWs + Dam + width_log + 
             no3_log + po4_log, # Adding nutrients in to be able to compare to
           # other covariates and previous model fit.
           random = ~1 | huc2_id, 
           method = "ML",
           data = dat_rmax_lm3)

plot(a9) # Looks alright.
summary(a9)

# Refit with REML for final full model

afinal_nuts <- lme(logrmax ~ summerT + NHD_RdDensWs + Dam + width_log + 
                     no3_log + po4_log,
              random = ~1 | huc2_id, 
              method = "REML",
              data = dat_rmax_lm3)

# Examine model diagnostics.
plot(afinal_nuts)
qqnorm(afinal_nuts) # Looking just fine.

# Examine model outputs.
aout_nuts <- summary(afinal_nuts)
coef(aout_nuts)

#                     Value  Std.Error DF    t-value      p-value
# (Intercept)  -1.379440083 0.35210413 67 -3.9177049 0.0002123089
# summerT       0.008003002 0.01281819 67  0.6243471 0.5345199231
# NHD_RdDensWs  0.012988349 0.01193486 67  1.0882695 0.2803753739
# Dam50        -0.167631532 0.15945979 67 -1.0512464 0.2969223232
# Dam80        -0.033620054 0.11455524 67 -0.2934833 0.7700598988
# Dam95         0.074057268 0.07251922 67  1.0212088 0.3108288147
# width_log     0.197495863 0.07636330 67  2.5862666 0.0118809167
# no3_log       0.163805166 0.08589546 67  1.9070295 0.0608054303
# po4_log       0.063028975 0.09046807 67  0.6966986 0.4884016515

##### Precipitation #####

# And build separate model for precipitation.

dat_rmax_lm4 <- dat_rmax_trim %>%
  dplyr::select(logrmax, summerT, width_log, NHD_RdDensWs, 
                Dam, huc2_id, pre_mm_cyr) %>%
  drop_na() # 115 sites left

a10 <- lme(logrmax ~ summerT + NHD_RdDensWs + Dam + width_log + 
            pre_mm_cyr, # Adding precip.
          random = ~1 | huc2_id, 
          method = "ML",
          data = dat_rmax_lm4)

plot(a10) # Looks alright.
qqnorm(a10) # Yep, totally fine.
summary(a10)

# Refit with REML for final full model

afinal_pre <- lme(logrmax ~ summerT + NHD_RdDensWs + Dam + width_log + 
                     pre_mm_cyr,
                   random = ~1 | huc2_id, 
                   method = "REML",
                   data = dat_rmax_lm4)

# Examine model diagnostics.
plot(afinal_pre)
qqnorm(afinal_pre) # Looking just fine.

# Examine model outputs.
aout_pre <- summary(afinal_pre)
coef(aout_pre)

#                      Value    Std.Error DF    t-value      p-value
# (Intercept)  -1.1196727973 0.3071445953 93 -3.6454257 0.0004396778
# summerT      -0.0161320206 0.0099431295 93 -1.6224289 0.1080968595
# NHD_RdDensWs  0.0379258201 0.0133298879 93  2.8451717 0.0054595347
# Dam50        -0.0520960691 0.1255499627 93 -0.4149429 0.6791393598
# Dam80        -0.0327519496 0.1047495372 93 -0.3126692 0.7552320717
# Dam95         0.0346755633 0.0730668281 93  0.4745733 0.6362036521
# width_log     0.3185124321 0.0835037112 93  3.8143506 0.0002455406
# pre_mm_cyr   -0.0000872829 0.0001621342 93 -0.5383375 0.5916302415

#### Model 2: Max. Algal Yield using 'brms' ####

# First, need to bind yield estimates with remaining data
# Using calculation based on eq. 7b from Scheuerell 2016 = `yield_med2`
# Combining with datasets from above.
dat_yield_combo <- left_join(dat_yield, site_precip, by = c("site_name" = "SiteID"))

dat_yield_combo <- left_join(dat_yield_combo, site_HUC2)

# And visualize the relationships with median yield values.
may_covs <- ggpairs(dat_yield_combo %>% 
                       dplyr::select(yield_med2, cvQ:Orthophosphate,
                                     pre_mm_cyr:huc2_id))

# ggsave(may_covs,
#        filename = "figures/teton_fall22/yield_covariates.jpg",
#        width = 50,
#        height = 50,
#        units = "cm")

# Some notes regarding these covariates.

# (1) I think I will again log transform the following variables: 
# width_med, Nitrate, and orthoP.

# (2) Will keep in mind correlations btw covariates found earlier.

# (3) The following jump out as potentially important (not including those
# data used to actually generate estimates): temperature, latitude, road
# density/pct impervious in watershed, empirical coefficient `a` of
# width estimation equation, maybe canals+dams??, width, precip, and highest
# outliers appear in certain HUC2s.

# Proposed starting model structure:

# max algal yield ~ temp + precip + roads + dams + width + 1 | HUC2

# Choosing precip over latitude for now, since they're highly
# correlated, but precip seems less skewed than the obvious geographic
# NE bias of these USGS gauges.

# Wil create a separate model adding in nutrients based on best model
# fit here (since records are far fewer).

# One on one plots for covariates of interest vs. may.

hist(dat_yield_combo$yield_med2)

# Going to log transform yield too.
dat_yield_combo <- dat_yield_combo %>%
  mutate(log_yield = log10(yield_med2)) %>%
  mutate(log_width = log10(width_med))

plot(log_yield ~ summerT, data = dat_yield_combo)
plot(log_yield ~ pre_mm_cyr, data = dat_yield_combo)
plot(log_yield ~ NHD_RdDensWs, data = dat_yield_combo)
plot(log_yield ~ Dam, data = dat_yield_combo)
plot(log_yield ~ log_width, data = dat_yield_combo)
plot(log_yield ~ huc2_id, data = dat_yield_combo)
hist(dat_yield_combo$log_yield)

# Ok, and making the final dataset with which to build models
dat_yield_brms <- dat_yield_combo %>%
  dplyr::select(log_yield, summerT, pre_mm_cyr,
                NHD_RdDensWs, Dam, log_width,
                huc2_id)

##### Step 1: Create multi-level model.

y1 <- brm(log_yield ~ summerT + pre_mm_cyr + NHD_RdDensWs +
            Dam + log_width + (1|huc2_id), 
         data = dat_yield_brms, family = gaussian())
# assumes 4 chains and 2000 iterations (1000 warm-up)
# started at 10:39am - finished at 10:40am :)

# Warning message:
# Rows containing NAs were excluded from the model. 

##### Step 2: Examine model outputs.

summary(y1)

#              Estimate Est.Error l-95% CI u-95% CI Rhat
# Intercept        0.55      0.37    -0.16     1.25 1.00
# summerT         -0.03      0.01    -0.05    -0.01 1.00
# pre_mm_cyr      -0.00      0.00    -0.00    -0.00 1.00
# NHD_RdDensWs     0.02      0.02    -0.01     0.05 1.00
# Dam50            0.10      0.16    -0.22     0.41 1.00
# Dam80            0.10      0.13    -0.16     0.36 1.00
# Dam95            0.19      0.09     0.02     0.37 1.00
# log_width        0.50      0.10     0.30     0.71 1.00

# Well, for one, this is great convergence! All Rhat < 1.05.
# And at first glance, temperature, dam95, and stream width jump out
# as important.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged well, so let's look at chain
# mixing and posterior distributions.
plot(y1, variable = c("b_summerT", "b_pre_mm_cyr", "b_NHD_RdDensWs",
                      "b_Dam50", "b_Dam80", "b_Dam95", "b_log_width"))

# Chains all appear well-mixed, but let's also check things in shinystan.
launch_shinystan(y1)

# No divergent transitions appear in the log posterior plots for any
# of the parameters.

# Finally, examine to be sure no n_eff are < 0.1
mcmc_plot(y1, type = "neff")

##### Step 4: Examine model relationships for each predictor.

plot(conditional_effects(y1, effects = "summerT"))
plot(conditional_effects(y1, effects = "pre_mm_cyr"))
plot(conditional_effects(y1, effects = "NHD_RdDensWs"))
plot(conditional_effects(y1, effects = "Dam"))
plot(conditional_effects(y1, effects = "log_width"))

# Note, can investigate scenarios like effect1:effect2 here
# and it will automatically choose percentiles to predict.

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
dat_yield_brms$obs <- c(1:length(dat_yield_brms$log_yield))

y1.1 <- brm(log_yield ~ summerT + pre_mm_cyr + NHD_RdDensWs +
            Dam + log_width + (1|huc2_id) + (1|obs), 
          data = dat_yield_brms, family = gaussian())
# 603 divergent transitions EEK!

# Compare with original model using leave-one-out approximation.
loo(y1, y1.1)

# Model comparisons:
#     elpd_diff se_diff
# y1.1   0.0       0.0  
# y1   -14.9       1.7 

# Higher expected log posterior density (elpd) values = better fit.

# So, in this case model accounting for overdispersion (y1.1) fits better.
# But there are 40 problematic observations...

# Compare with WAIC instead.
waic1 <- waic(y1)
waic2 <- waic(y1.1)
loo_compare(waic1, waic2)

# Same output as above.
#      elpd_diff se_diff
# y1.1   0.0       0.0  
# y1   -17.4       1.2

##### Step 6: Plot the results.

get_variables(y1)

# b_Intercept refers to global mean
# r_huc2_id[] are the offsets from that mean for each condition

(y_fig <- mcmc_plot(y1, variable = c("b_summerT", "b_pre_mm_cyr",
                                    "b_NHD_RdDensWs", "b_Dam50", 
                                    "b_Dam80", "b_Dam95", "b_log_width"),
      #type = "intervals",
      point_est = "median", # default = "median"
      prob = 0.66, # default = 0.5
      prob_outer = 0.95) + # default = 0.9
    vline_at(v = 0) +
    labs(x = "Posterior",
         y = "Coefficients") +
    theme_bw())

# Save out this figure.
# ggsave(y_fig,
#        filename = "figures/teton_fall22/brms_yield_021423.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

# Can also use pars = c("^r_", "^b_", "^sd_") in place of variable phrasing
# to see all results.

# Default is type = "intervals".

# Also trying to figure out a way to plot histograms/density plots
# over top of these intervals.

my_posterior <- as.matrix(y1)

(y_fig2 <- mcmc_areas(my_posterior,
           # removed precip bc it was flattening everything
           # then removed all but "significant" covariates
           pars = c("b_summerT", "b_Dam95", "b_log_width"),
           prob = 0.95) +
  ggtitle("Posterior distributions",
          "with 95% credible intervals not crossing zero") +
  vline_at(v = 0) +
  labs(x = "Posterior",
       y = "Coefficients") +
  theme_bw())

# Save out this figure.
# ggsave(y_fig2,
#        filename = "figures/teton_fall22/brms_yield_sig_021423.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

##### Latitude #####

# Proposed starting model structure:

# max algal yield ~ temp + lat + roads + dams + width + 1 | HUC2

# Choosing latitude here, since it drops less observations.

# One on one plots for covariates of interest vs. may.

hist(dat_yield_combo$yield_med2)
plot(log_yield ~ Lat_WGS84, data = dat_yield_combo)

# Ok, and making the final dataset with which to build models
dat_yield_brms2 <- dat_yield_combo %>%
  dplyr::select(log_yield, summerT, Lat_WGS84,
                NHD_RdDensWs, Dam, log_width,
                huc2_id)

##### Step 1: Create multi-level model.

y2 <- brm(log_yield ~ summerT + Lat_WGS84 + NHD_RdDensWs +
            Dam + log_width + (1|huc2_id), 
          data = dat_yield_brms2, family = gaussian())
# assumes 4 chains and 2000 iterations (1000 warm-up)
# started at 10:53am - finished at 10:54am :)

# Warning message:
# Rows containing NAs were excluded from the model. 

##### Step 2: Examine model outputs.

summary(y2)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept        0.38      0.84    -1.29     2.05 1.00     3285     2745
# summerT         -0.03      0.02    -0.06    -0.00 1.00     4034     3519
# Lat_WGS84       -0.00      0.02    -0.03     0.03 1.00     3104     2767
# NHD_RdDensWs     0.02      0.02    -0.01     0.05 1.00     3543     2774
# Dam50           -0.04      0.13    -0.30     0.21 1.00     5807     3232
# Dam80            0.17      0.11    -0.05     0.40 1.00     5305     2545
# Dam95            0.18      0.08     0.02     0.34 1.00     4559     3233
# log_width        0.50      0.10     0.31     0.69 1.00     3431     3238

# Well, great convergence again! All Rhat < 1.05.
# And temperature, dam95, and stream width again jump out.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged well, so let's look at chain
# mixing and posterior distributions.
plot(y2, variable = c("b_summerT", "b_Lat_WGS84", "b_NHD_RdDensWs",
                      "b_Dam50", "b_Dam80", "b_Dam95", "b_log_width"))

# Chains all appear well-mixed, but let's also check things in shinystan.
launch_shinystan(y2)

# No divergent transitions appear in the log posterior plots for any
# of the parameters.

# Finally, examine to be sure no n_eff are < 0.1
mcmc_plot(y2, type = "neff")

##### Step 4: Examine model relationships for each predictor.

plot(conditional_effects(y2, effects = "summerT"))
plot(conditional_effects(y2, effects = "Lat_WGS84"))
plot(conditional_effects(y2, effects = "NHD_RdDensWs"))
plot(conditional_effects(y2, effects = "Dam"))
plot(conditional_effects(y2, effects = "log_width"))

# Note, can investigate scenarios like effect1:effect2 here
# and it will automatically choose percentiles to predict.

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
dat_yield_brms2$obs <- c(1:length(dat_yield_brms2$log_yield))

y2.1 <- brm(log_yield ~ summerT + Lat_WGS84 + NHD_RdDensWs +
              Dam + log_width + (1|huc2_id) + (1|obs), 
            data = dat_yield_brms2, family = gaussian())
# 76 divergent transitions

# Compare with original model using leave-one-out approximation.
loo(y2, y2.1)

# Model comparisons:
#     elpd_diff se_diff
# y2.1   0.0       0.0  
# y2    -4.4       0.7 

# Higher expected log posterior density (elpd) values = better fit.

# So, in this case model accounting for overdispersion (y2.1) fits better.
# But there are 22 problematic observations...

##### Step 6: Plot the results.

get_variables(y2)

# b_Intercept refers to global mean
# r_huc2_id[] are the offsets from that mean for each condition

(y2_fig <- mcmc_plot(y2, variable = c("b_summerT", "b_Lat_WGS84",
                                     "b_NHD_RdDensWs", "b_Dam50", 
                                     "b_Dam80", "b_Dam95", "b_log_width"),
                    #type = "intervals",
                    point_est = "median", # default = "median"
                    prob = 0.66, # default = 0.5
                    prob_outer = 0.95) + # default = 0.9
    vline_at(v = 0) +
    labs(x = "Posterior",
         y = "Coefficients") +
    theme_bw())

# Save out this figure.
# ggsave(y2_fig,
#        filename = "figures/teton_fall22/brms_yield2_021623.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

# Also trying to figure out a way to plot histograms/density plots
# over top of these intervals.

my_posterior2 <- as.matrix(y2)

(y2_fig2 <- mcmc_areas(my_posterior2,
                      # removed lat bc it was flattening everything
                      # then removed all but "significant" covariates
                      pars = c("b_summerT", "b_Dam95", "b_log_width"),
                      prob = 0.95) +
    ggtitle("Posterior distributions",
            "with 95% credible intervals not crossing zero") +
    vline_at(v = 0) +
    labs(x = "Posterior",
         y = "Coefficients") +
    theme_bw())

# Save out this figure.
# ggsave(y2_fig2,
#        filename = "figures/teton_fall22/brms_yield2_sig_021623.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

##### Nutrients #####

# And build separate model for nutrients.
# Proposed starting model structure:

# max algal yield ~ temp + lat + roads + dams + width + NO3 + PO4 + 1 | HUC2

# Sticking with latitude here, since it drops less observations.

# Need to log transform nutrients.
dat_yield_combo <- dat_yield_combo %>%
  mutate(no3_log = log10(Nitrate),
         po4_log = log10(Orthophosphate))

# One on one plots for covariates of interest vs. may.

plot(log_yield ~ no3_log, data = dat_yield_combo)
plot(log_yield ~ po4_log, data = dat_yield_combo)

# Ok, and making the final dataset with which to build models
dat_yield_brms3 <- dat_yield_combo %>%
  dplyr::select(log_yield, summerT, Lat_WGS84,
                NHD_RdDensWs, Dam, log_width,
                no3_log, po4_log, huc2_id)

##### Step 1: Create multi-level model.

y3 <- brm(log_yield ~ summerT + Lat_WGS84 + NHD_RdDensWs +
            Dam + log_width + no3_log + po4_log + (1|huc2_id), 
          data = dat_yield_brms3, family = gaussian(),
          iter = 3000) 
# upped to 3k, bc 2 div. transitions with 2000 iterations
# started at 11:17am - finished at 11:20am :)

# 1 divergent transition.

##### Step 2: Examine model outputs.

summary(y3)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept        0.40      1.24    -2.08     2.78 1.00     4411     3891
# summerT         -0.02      0.03    -0.07     0.03 1.00     4646     4515
# Lat_WGS84        0.00      0.02    -0.04     0.04 1.00     4243     3839
# NHD_RdDensWs    -0.01      0.02    -0.04     0.03 1.00     5928     4530
# Dam50           -0.20      0.21    -0.61     0.21 1.00     8027     4546
# Dam80            0.01      0.15    -0.29     0.31 1.00     9333     4600
# Dam95            0.13      0.10    -0.06     0.32 1.00     8336     4952
# log_width        0.37      0.12     0.14     0.60 1.00     4768     4785
# no3_log          0.02      0.12    -0.21     0.25 1.00     5186     4633
# po4_log          0.14      0.12    -0.10     0.38 1.00     5485     4415

# Well, great convergence again! All Rhat < 1.05.
# And stream width again jumps out.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged well, so let's look at chain
# mixing and posterior distributions.
plot(y3, variable = c("b_summerT", "b_Lat_WGS84", "b_NHD_RdDensWs",
                      "b_Dam50", "b_Dam80", "b_Dam95", "b_log_width",
                      "b_no3_log", "b_po4_log"))

# Chains all appear well-mixed, but let's also check things in shinystan.
launch_shinystan(y3)

# 1 divergent transition.

# Finally, examine to be sure no n_eff are < 0.1
mcmc_plot(y3, type = "neff")

##### Step 4: Examine model relationships for each predictor.

plot(conditional_effects(y3, effects = "summerT"))
plot(conditional_effects(y3, effects = "Lat_WGS84"))
plot(conditional_effects(y3, effects = "NHD_RdDensWs"))
plot(conditional_effects(y3, effects = "Dam"))
plot(conditional_effects(y3, effects = "log_width"))
plot(conditional_effects(y3, effects = "no3_log"))
plot(conditional_effects(y3, effects = "po4_log"))

# Note, can investigate scenarios like effect1:effect2 here
# and it will automatically choose percentiles to predict.

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
dat_yield_brms3$obs <- c(1:length(dat_yield_brms3$log_yield))

y3.1 <- brm(log_yield ~ summerT + Lat_WGS84 + NHD_RdDensWs +
              Dam + log_width + no3_log + po4_log +
              (1|huc2_id) + (1|obs), 
            data = dat_yield_brms3, family = gaussian())
# 13 divergent transitions

# Compare with original model using leave-one-out approximation.
loo(y3, y3.1)

# Model comparisons:
#     elpd_diff se_diff
# y3.1   0.0       0.0  
# y3   -11.3       1.0 

# Higher expected log posterior density (elpd) values = better fit.

# So, in this case model accounting for overdispersion (y3.1) fits better.
# But there are 33 problematic observations...

##### Step 6: Plot the results.

get_variables(y3)

# b_Intercept refers to global mean
# r_huc2_id[] are the offsets from that mean for each condition

(y3_fig <- mcmc_plot(y3, variable = c("b_summerT", "b_Lat_WGS84",
                                      "b_NHD_RdDensWs", "b_Dam50", 
                                      "b_Dam80", "b_Dam95", "b_log_width",
                                      "b_no3_log", "b_po4_log"),
                     #type = "intervals",
                     point_est = "median", # default = "median"
                     prob = 0.66, # default = 0.5
                     prob_outer = 0.95) + # default = 0.9
    vline_at(v = 0) +
    labs(x = "Posterior",
         y = "Coefficients") +
    theme_bw())

# Save out this figure.
# ggsave(y3_fig,
#        filename = "figures/teton_fall22/brms_yield3_021623.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

# Also trying to figure out a way to plot histograms/density plots
# over top of these intervals.

my_posterior3 <- as.matrix(y3)

(y3_fig2 <- mcmc_areas(my_posterior3,
                       # removed lat bc it was flattening everything
                       # then removed all but "significant" covariates
                       pars = c("b_log_width"),
                       prob = 0.95) +
    ggtitle("Posterior distributions",
            "with 95% credible intervals not crossing zero") +
    vline_at(v = 0) +
    labs(x = "Posterior",
         y = "Coefficients") +
    theme_bw())

# Save out this figure.
# ggsave(y3_fig2,
#        filename = "figures/teton_fall22/brms_yield3_sig_021623.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

#### Model 3: Qc:Q2yr ####

# Trim imported data down to variables of interest.

dat_Qc_trim <- dat_Qc %>%
  select(site_name,
         c_med, Qc_Q2yr, cvQ, meanGPP:NHD_AREASQKM, 
         NHD_RdDensCat:NHD_PctImp2011Ws, 
         Canal:width_med)

# Join with precip data.
dat_Qc_trim <- left_join(dat_Qc_trim, site_precip, by = c("site_name" = "SiteID"))

# Join with HUC2 data.
dat_Qc_trim <- left_join(dat_Qc_trim, site_HUC2)

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
# than they did in the rmax covariate exploration. So, makes double sense
# to remove them.

# (3) Remaining Pearson's correlation values are below 0.5.

# Log transform necessary covariates.

dat_Qc_trim <- dat_Qc_trim %>%
  mutate(GPP_log = log10(meanGPP),
         area_log = log10(NHD_AREASQKM),
         width_log = log10(width_med))

# Notes on model structure:

# Qc:Q2yr ~ precip + size + roads + dams + 1 | HUC2

# I am going to include the following covariates as representatives of the
# corresponding environmental factors:

# precipitation - land and water connectivity metric
# order - stream size**
# watershed area - stream size**
# width - stream size**
# road density - terrestrial development
# dam - aquatic development
# HUC2 - to account for random effect of geography

# ** Models will explore each of the size indicators separately. Although
# I lean towards using width to mirror rmax above.

# The following covariates have been removed for the following reasons:

# light - not a factor for flow disturbance thresholds
# temperature - not a factor for flow disturbance thresholds
# cvQ - flow but used to estimate Qc
# GPP - biological productivity (function of flow & light) but used to 
# estimate Qc
# latitude - a rough location/temperature index?
# longitude - a rough location/aridity index?
# % impervious land cover - too closely correlated with road density
# canal - another metric of terr/aq development but felt duplicative
# NO3/PO4 - not a factor for flow disturbance thresholds

# One on one plots for covariates of interest vs. QcQ2yr.
hist(dat_Qc_trim$Qc_Q2yr)

# Going to log transform QcQ2yr too.
dat_Qc_trim <- dat_Qc_trim %>%
  mutate(logQcQ2max = log10(Qc_Q2yr))

plot(logQcQ2max ~ area_log, data = dat_Qc_trim)
plot(logQcQ2max ~ Order, data = dat_Qc_trim) # Look about the same
plot(logQcQ2max ~ width_log, data = dat_Qc_trim)
plot(logQcQ2max ~ NHD_RdDensWs, data = dat_Qc_trim)
plot(logQcQ2max ~ Dam, data = dat_Qc_trim)
plot(logQcQ2max ~ pre_mm_cyr, data = dat_Qc_trim)
plot(logQcQ2max ~ huc2_id, data = dat_Qc_trim)
hist(dat_Qc_trim$logQcQ2max)

# Ok, and making the final dataset with which to build models since I can't have
# NAs with many of these functions.
dat_Qc_lm <- dat_Qc_trim %>%
  dplyr::select(logQcQ2max, Order, area_log, 
                width_log, NHD_RdDensWs, Dam, huc2_id) %>%
  drop_na() # 151 sites left

##### Step 1: Create lm() and check residuals.

# rmax ~ temp + roads + dams

c1 <- lm(logQcQ2max ~ NHD_RdDensWs + Dam, 
         data = dat_Qc_lm)

plot(c1) # looks alright

summary(c1) # examine initial model output without the grouping by watershed
# hmmm, let's see how this does with the other covariates added in

##### Step 2: Fit the lm() with GLS and compare to lme().

c2 <- gls(logQcQ2max ~ NHD_RdDensWs + Dam, 
          data = dat_Qc_lm) # effectively a lm

c3 <- lme(logQcQ2max ~ NHD_RdDensWs + Dam, 
          random = ~1 | huc2_id, data = dat_Qc_lm) # with random effect

anova(c2, c3) # c2 preferred, but I want to keep the random term in

##### Step 3: Decide on a variance structure.

# None at this time.

##### Step 4,5,6: Fit the lme(), compare with lm(), and check residuals.

# See Steps 2 & 3.

##### Step 7/8: Step-wise Optimal Fixed Structure.

# Try different covariates for stream size.
c4 <- lme(logQcQ2max ~ NHD_RdDensWs + Dam + Order, 
          random = ~1 | huc2_id, 
          method = "ML",
          data = dat_Qc_lm) # stream order

c5 <- lme(logQcQ2max ~ NHD_RdDensWs + Dam + area_log, 
          random = ~1 | huc2_id, 
          method = "ML",
          data = dat_Qc_lm) # watershed area

c6 <- lme(logQcQ2max ~ NHD_RdDensWs + Dam + width_log, 
          random = ~1 | huc2_id, 
          method = "ML",
          data = dat_Qc_lm) # stream width

anova(c4, c5, c6) # order (c5) AIC is least, but moving forward with width
# since it's better than order (which is somewhat subjective).

# Investigate precipitation, but keep in mind, it drops 40+ records.
# (Archiving this because investigating this below.)
# c6.2 <- lme(logQcQ2max ~ NHD_RdDensWs + Dam + width_log, 
#             random = ~1 | huc2_id, 
#             method = "ML",
#             data = dat_Qc_trim %>%
#               dplyr::select(logQcQ2max, NHD_RdDensWs, Dam,
#                             width_log, pre_mm_cyr, huc2_id) %>%
#               drop_na())

# c7 <- lme(logQcQ2max ~ NHD_RdDensWs + Dam + width_log +
#             pre_mm_cyr, 
#           random = ~1 | huc2_id, 
#           method = "ML",
#           data = dat_Qc_trim %>%
#             dplyr::select(logQcQ2max, NHD_RdDensWs, Dam,
#                           width_log, pre_mm_cyr, huc2_id) %>%
#             drop_na())

#anova(c6.2, c7) # Identical, so since it forces us to drop so many records 
# (n=115 remaining), I'm skipping it.

##### Step 9: Refit with REML for final full model

# Maximize data availability:
dat_Qc_lm2 <- dat_Qc_trim %>%
  dplyr::select(logQcQ2max, width_log, NHD_RdDensWs, 
                Dam, huc2_id) %>%
  drop_na() #151

# Precip didn't seem to add anything to the model, and subtracted >40 records,
# so choosing not to use it at this time.
cfinal <- lme(logQcQ2max ~ # log-transformed because right-skewed
                NHD_RdDensWs + # used by JB as metric of development, more sig. than Cat
                Dam + # more pertinent to our interests than "Canal"
                width_log, # not best performing, but more trustworthy than order
              random = ~1 | huc2_id, # because HUC10 was too few sites/watershed
              method = "REML",
              data = dat_Qc_lm2) # n = 133

# Examine model diagnostics.
plot(cfinal) # YAY! Residuals looking great.
qqnorm(cfinal) # And qqplot looks just fine.

# Examine model outputs.
cout <- summary(cfinal)
coef(cout)

#                   Value  Std.Error  DF    t-value   p-value
# (Intercept)  -0.23264861 0.19330081 113 -1.2035573 0.2312760
# NHD_RdDensWs -0.02130914 0.01818270 113 -1.1719462 0.2436836
# Dam50         0.29308797 0.18682568 113  1.5687777 0.1194958
# Dam80         0.05513241 0.14634429 113  0.3767309 0.7070801
# Dam95         0.02196567 0.09679328 113  0.2269339 0.8208852
# width_log    -0.01459310 0.10789160 113 -0.1352570 0.8926493

r.squaredGLMM(cfinal)

#             R2m        R2c
# [1,] 0.02458802 0.06703832

##### Precipitation #####

# Maximize data availability:
dat_Qc_lm3 <- dat_Qc_trim %>%
  dplyr::select(logQcQ2max, width_log, NHD_RdDensWs, 
                Dam, huc2_id, pre_mm_cyr) %>%
  drop_na() #103

c8 <- lme(logQcQ2max ~ NHD_RdDensWs + Dam + width_log + 
             pre_mm_cyr, # Adding precip.
           random = ~1 | huc2_id, 
           method = "ML",
           data = dat_Qc_lm3)

plot(c8) # Looks alright.
qqnorm(c8) # Yep, totally fine.
summary(c8)

# Refit with REML for final full model

cfinal_pre <- lme(logQcQ2max ~ NHD_RdDensWs + Dam + width_log + 
                    pre_mm_cyr,
                  random = ~1 | huc2_id, 
                  method = "REML",
                  data = dat_Qc_lm3)

# Examine model diagnostics.
plot(cfinal_pre)
qqnorm(cfinal_pre) # Looking just fine.

# Examine model outputs.
cout_pre <- summary(cfinal_pre)
coef(cout_pre)

#                      Value    Std.Error DF      t-value   p-value
# (Intercept)   0.0930882014 0.3122302087 82  0.298139638 0.7663507
# NHD_RdDensWs -0.0126805666 0.0194349217 82 -0.652462963 0.5159275
# Dam50         0.3012597661 0.2146412775 82  1.403550005 0.1642286
# Dam80        -0.0017105104 0.1756896271 82 -0.009735978 0.9922556
# Dam95        -0.0838730873 0.1093675140 82 -0.766892145 0.4453480
# width_log    -0.0592449483 0.1150885535 82 -0.514777069 0.6080936
# pre_mm_cyr   -0.0002641592 0.0002098016 82 -1.259090714 0.2115701

# End of script.