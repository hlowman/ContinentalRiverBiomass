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
lapply(c("tidybayes", "brms", "tidyverse", "lubridate", 
         "data.table", "GGally",
         "multcomp", "patchwork", "bayesplot",
         "modelsummary", "here", "nlme","loo"), 
       require, character.only=T)

#### Data ####

# Import necessary datasets.

# First, the data for the rmax models.
dat_rmax <- readRDS("data_working/rmax_filtered_159sites_113022.rds")

# Next, the data for the Qc:Q2yr models.
dat_Qc <- readRDS("data_working/QcQ2_filtered_141sites_113022.rds")

# Also , the data for maximum algal yields.
dat_yield <- readRDS("data_working/maxalgalyield_159sites_021323.rds")

# Finally, the data for sites' Qc exceedances.
dat_exc <- readRDS("data_working/Qc_exceedances_159sites_021423.rds")

# And the hypoxia dataset for additional info re: precip. "pre_mm_cyr"
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")

# And the dataset with all HUC delineations.
site_HUC <- readRDS("data_working/HUC12_159sites_120922.rds")

# And the newer nutrient dataset.
new_nut <- readRDS("data_working/USGS_WQP_nuts_aggsite_022322.rds")

# Select for additional variables of interest.

# Annual average precip for local catchment/watershed from HydroATLAS (mm)
site_precip <- site_info %>%
  dplyr::select(SiteID, pre_mm_cyr, pre_mm_uyr)

site_HUC2 <- site_HUC %>%
  dplyr::select(site_name, huc2_id)

site_P <- new_nut %>%
  filter(CharacteristicName == "Phosphorus") %>%
  rename(Phosphorus = "mean_mg_L")

site_P$site_name <- str_replace_all(site_P$MonitoringLocationIdentifier, 'USGS-', 'nwis_')

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
dat_yield_combo <- left_join(dat_yield, site_precip, 
                             by = c("site_name" = "SiteID"))

dat_yield_combo <- left_join(dat_yield_combo, site_HUC2)

dat_yield_combo <- left_join(dat_yield_combo, site_P)

# And visualize the relationships with median yield values.
may_covs <- ggpairs(dat_yield_combo %>% 
                       dplyr::select(yield_med2, cvQ:Orthophosphate,
                                     Phosphorus,
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

# max algal yield/accrual ~ size + roads + dams + temperature + 1 | HUC2

# Removed latitude/precip following conversations with Joanna.

# Wil create a separate model adding in nutrients based on best model
# fit here (since records are far fewer).

# One on one plots for covariates of interest vs. may.

hist(dat_yield_combo$yield_med2)

# Going to log transform yield too.
dat_yield_combo <- dat_yield_combo %>%
  mutate(log_yield = log10(yield_med2)) %>%
  mutate(log_width = log10(width_med)) %>%
  # Also creating a new categorical dam column to model by.
  # same metric used in QC:Q2yr model below.
  mutate(Dam_binary = factor(case_when(
    Dam %in% c("50", "80", "95") ~ "0", # Potential = 5-50%
    Dam == "0" ~ "1", # Certain = 100%
    TRUE ~ NA)))

plot(log_yield ~ summerT, data = dat_yield_combo)
plot(log_yield ~ NHD_RdDensWs, data = dat_yield_combo)
plot(log_yield ~ Dam_binary, data = dat_yield_combo)
plot(log_yield ~ log_width, data = dat_yield_combo)
plot(log_yield ~ huc2_id, data = dat_yield_combo)
hist(dat_yield_combo$log_yield)

# Ok, and making the final dataset with which to build models
dat_yield_brms <- dat_yield_combo %>%
  dplyr::select(log_yield, summerT, NHD_RdDensWs, 
                Dam_binary, log_width, huc2_id)

##### Step 1: Create multi-level model.

y1 <- brm(log_yield ~ log_width + NHD_RdDensWs +
            Dam_binary + summerT + (1|huc2_id), 
         data = dat_yield_brms, family = gaussian())
# assumes 4 chains and 2000 iterations (1000 warm-up)
# started at 10:39am - finished at 10:40am on the server :)

# Export for safekeeping.
# saveRDS(y1, "data_posthoc_modelfits/accrual_brms_030123.rds")

##### Step 2: Examine model outputs.

summary(y1)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept        0.51      0.29    -0.05     1.08 1.00     4481     2960
# log_width        0.46      0.09     0.29     0.65 1.00     3972     3159
# NHD_RdDensWs     0.02      0.02    -0.01     0.05 1.00     4094     3222
# Dam_binary1     -0.14      0.07    -0.29     0.00 1.00     5189     3097
# summerT         -0.03      0.01    -0.06    -0.01 1.00     4332     2907

# Well, for one, this is great convergence! All Rhat < 1.05.
# And at first glance, temperature, dam_binary1, and stream width jump out
# as important.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged well, so let's look at chain
# mixing and posterior distributions.
plot(y1, variable = c("b_summerT", "b_NHD_RdDensWs",
                      "b_Dam_binary1", "b_log_width"))

# Chains all appear well-mixed, but let's also check things in shinystan.
launch_shinystan(y1)

# No divergent transitions appear in the log posterior plots for any
# of the parameters.

# Finally, examine to be sure no n_eff are < 0.1
mcmc_plot(y1, type = "neff")

##### Step 4: Examine model relationships for each predictor.

plot(conditional_effects(y1, effects = "log_width"))
plot(conditional_effects(y1, effects = "NHD_RdDensWs"))
plot(conditional_effects(y1, effects = "Dam_binary"))
plot(conditional_effects(y1, effects = "summerT"))

# Note, can investigate scenarios like effect1:effect2 here
# and it will automatically choose percentiles to predict.

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
dat_yield_brms$obs <- c(1:length(dat_yield_brms$log_yield))

y1.1 <- brm(log_yield ~ log_width + NHD_RdDensWs +
            Dam_binary + summerT + (1|huc2_id) + (1|obs), 
          data = dat_yield_brms, family = gaussian())
# 316 divergent transitions EEK!

# Compare with original model using leave-one-out approximation.
loo(y1, y1.1)

# Model comparisons:
#     elpd_diff se_diff
# y1.1   0.0       0.0  
# y1   -10.4       1.0 

# Higher expected log posterior density (elpd) values = better fit.
# So, in this case model accounting for overdispersion (y1.1) fits better.
# But there are 30 problematic observations...
# And a few hundred divergences and poor Rhat values, so my gut says the
# original model (y1) has the better fit.

##### Step 6: Plot the results.

get_variables(y1)

# b_Intercept refers to global mean
# r_huc2_id[] are the offsets from that mean for each condition

(y_fig <- mcmc_plot(y1, variable = c("b_log_width", "b_NHD_RdDensWs",
                                     "b_Dam_binary1", "b_summerT"),
      #type = "intervals",
      point_est = "median", # default = "median"
      prob = 0.66, # default = 0.5
      prob_outer = 0.95) + # default = 0.9
    vline_at(v = 0) +
    labs(x = "Posterior Estimates",
         y = "Predictors") +
    scale_y_discrete(labels = c("b_log_width" = "log(Width)",
                                "b_NHD_RdDensWs" = "Road Density",
                                "b_Dam_binary1" = "Dam",
                                "b_summerT" = "Temperature")) +
    theme_bw() +
    theme(text = element_text(size = 16)))

# Save out this figure.
# ggsave(y_fig,
#        filename = "figures/teton_fall22/brms_yield_030123.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

# Can also use pars = c("^r_", "^b_", "^sd_") in place of variable phrasing
# to see all results.

# Default is type = "intervals".

# Also create conditional plots for each parameter.
# Plot conditional effects of all covariates.
# Using code from here to make them ggplots:
# https://bookdown.org/content/4857/conditional-manatees.html#summary-bonus-conditional_effects

y_t <- conditional_effects(y1, effects = "summerT")

# Create new dataframe
yt_df <- y_t$summerT %>%
  # and calculate true yield values
  # checked to be sure these default probs are 0.95
  mutate(yield = 10^`estimate__`,
         loweryield = 10^`lower__`,
         upperyield = 10^`upper__`)

(plot_yt <- ggplot(yt_df, aes(x = summerT, y = yield)) +
    geom_line(color = "black", linewidth = 1) +
    geom_ribbon(aes(ymin = loweryield, ymax = upperyield),
                alpha = 0.25) +
    geom_point(data = dat_yield_brms, aes(x = summerT, y = 10^log_yield),
               alpha = 0.2) +
    labs(x = paste0("Mean Summer Temperature (", '\u00B0', "C)"),
         y = expression(a[max])) +
    ylim(0, 17) +
    theme_bw())

y_d <- conditional_effects(y1, effects = "Dam_binary")

# Create new dataframe
yd_df <- y_d$Dam_binary %>%
  # and calculate true yield values
  mutate(yield = 10^`estimate__`,
         loweryield = 10^`lower__`,
         upperyield = 10^`upper__`)

(plot_yd <- ggplot(yd_df, aes(x = Dam_binary, y = yield)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = loweryield, ymax = upperyield), width = 0.2) +
    geom_jitter(data = dat_yield_brms, aes(x = Dam_binary, y = 10^log_yield),
               alpha = 0.1, width = 0.1) +
    labs(x = "Likelihood of Interference by Dams",
         y = expression(a[max])) +
    scale_x_discrete(labels = c("5-50%", "100%")) +
    ylim(0, 17) +
    theme_bw())

y_w <- conditional_effects(y1, effects = "log_width")

# Create new dataframe
yw_df <- y_w$log_width %>%
  # and calculate true yield values
  mutate(yield = 10^`estimate__`,
         loweryield = 10^`lower__`,
         upperyield = 10^`upper__`)

(plot_yw <- ggplot(yw_df, aes(x = 10^log_width, y = yield)) +
    geom_line(color = "black", linewidth = 1) +
    geom_ribbon(aes(ymin = loweryield, ymax = upperyield),
                alpha = 0.25) +
    geom_point(data = dat_yield_brms, 
               aes(x = 10^log_width, y = 10^log_yield),
               alpha = 0.2) +
    scale_x_log10()+
    labs(x = "River Width (m)",
         y = expression(a[max])) +
    ylim(0, 17) +
    theme_bw())

y_rd <- conditional_effects(y1, effects = "NHD_RdDensWs")

# Create new dataframe
yrd_df <- y_rd$NHD_RdDensWs %>%
  # and calculate true yield values
  mutate(yield = 10^`estimate__`,
         loweryield = 10^`lower__`,
         upperyield = 10^`upper__`)

(plot_yrd <- ggplot(yrd_df, aes(x = NHD_RdDensWs, y = yield)) +
    geom_line(color = "black", linewidth = 1) +
    geom_ribbon(aes(ymin = loweryield, ymax = upperyield),
                alpha = 0.25) +
    geom_point(data = dat_yield_brms, 
               aes(x = NHD_RdDensWs, y = 10^log_yield),
               alpha = 0.2) +
    #scale_x_log10()+
    labs(x = expression(Road~Density~by~Watershed~(km/km^2)),
         y = expression(a[max])) +
    ylim(0, 17) +
    theme_bw())

# Now, let's combine the above using patchwork.
(fig_cond_yield <- (y_fig | (plot_yt + plot_yd) / (plot_yrd + plot_yw)) +
    plot_annotation(tag_levels = 'A'))

# And export.
# ggsave(fig_cond_yield,
#        filename = "figures/teton_fall22/brms_yield_cond_032323.jpg",
#        width = 40,
#        height = 20,
#        units = "cm")

##### Nutrients #####

# And build separate model for nutrients.
# Proposed starting model structure:

# max algal yield ~ temp + roads + dams + width + NO3 + P + 1 | HUC2

# Need to log transform nutrients.
dat_yield_combo <- dat_yield_combo %>%
  mutate(no3_log = log10(Nitrate),
         p_log = log10(Phosphorus))

# One on one plots for covariates of interest vs. may.
plot(log_yield ~ no3_log, data = dat_yield_combo)
plot(log_yield ~ p_log, data = dat_yield_combo)

# Ok, and making the final dataset with which to build models
dat_yield_brms3 <- dat_yield_combo %>%
  dplyr::select(log_yield, summerT, NHD_RdDensWs,
                Dam_binary, log_width, no3_log, p_log, 
                huc2_id)

# And need to drop NAs.
dat_yield_brms3 <- dat_yield_brms3 %>%
  drop_na(no3_log, p_log) # 60 observations *sigh*

##### Step 1: Create multi-level model.

y3 <- brm(log_yield ~ summerT + NHD_RdDensWs + Dam_binary +
            log_width + no3_log + p_log + (1|huc2_id), 
          data = dat_yield_brms3, family = gaussian(),
          iter = 2000) # 1 divergent transition

# Export model fit for safekeeping.
# saveRDS(y3, "data_posthoc_modelfits/accrual_nuts_brms_030123.rds")

##### Step 2: Examine model outputs.

summary(y3)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept       -0.15      0.60    -1.40     1.01 1.00     2421     1320
# summerT          0.01      0.02    -0.04     0.06 1.00     2429     1572
# NHD_RdDensWs     0.00      0.02    -0.04     0.05 1.00     1930     1548
# Dam_binary1     -0.12      0.11    -0.33     0.09 1.00     3445     2615
# log_width        0.43      0.13     0.18     0.70 1.00     2296     1844
# no3_log         -0.05      0.14    -0.32     0.23 1.00     2831     2309
# p_log            0.19      0.15    -0.10     0.48 1.00     2417     2588

# Well, great convergence again! All Rhat < 1.05.
# And stream width again jumps out. But everything else seems to have
# been swamped out by the addition of nutrients.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged well, so let's look at chain
# mixing and posterior distributions.
plot(y3, variable = c("b_summerT", "b_NHD_RdDensWs",
                      "b_Dam_binary1", "b_log_width",
                      "b_no3_log", "b_p_log"))

# Chains all appear well-mixed, but let's also check things in shinystan.
launch_shinystan(y3)
# 1 divergent transition.

# Finally, examine to be sure no n_eff are < 0.1
mcmc_plot(y3, type = "neff")

##### Step 4: Examine model relationships for each predictor.

plot(conditional_effects(y3, effects = "summerT"))
plot(conditional_effects(y3, effects = "NHD_RdDensWs"))
plot(conditional_effects(y3, effects = "Dam_binary"))
plot(conditional_effects(y3, effects = "log_width"))
plot(conditional_effects(y3, effects = "no3_log"))
plot(conditional_effects(y3, effects = "p_log"))

# Note, can investigate scenarios like effect1:effect2 here
# and it will automatically choose percentiles to predict.

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
dat_yield_brms3$obs <- c(1:length(dat_yield_brms3$log_yield))

y3.1 <- brm(log_yield ~ summerT + NHD_RdDensWs +
              Dam_binary + log_width + no3_log + p_log +
              (1|huc2_id) + (1|obs), 
            data = dat_yield_brms3, family = gaussian())
# 146 divergent transitions eeeegads

# Compare with original model using leave-one-out approximation.
loo(y3, y3.1)

# Model comparisons:
#     elpd_diff se_diff
# y3.1   0.0       0.0  
# y3   -15.2       1.2 

# Higher expected log posterior density (elpd) values = better fit.

# So, in this case model accounting for overdispersion (y3.1) fits better.
# But there are 26 problematic observations, and using the logic I
# employed above, I'm sticking with the original model since it
# had far fewer divergent transitions.

##### Step 6: Plot the results.

get_variables(y3)

# b_Intercept refers to global mean
# r_huc2_id[] are the offsets from that mean for each condition

(y3_fig <- mcmc_plot(y3, variable = c("b_log_width", "b_p_log",
                                      "b_summerT", "b_NHD_RdDensWs",
                                      "b_no3_log", "b_Dam_binary1"),
                     #type = "intervals",
                     point_est = "median", # default = "median"
                     prob = 0.66, # default = 0.5
                     prob_outer = 0.95) + # default = 0.9
    vline_at(v = 0) +
    labs(x = "Posterior Estimates",
         y = "Predictors") +
    scale_y_discrete(labels = c("b_p_log" = "Phosphorus",
                                "b_no3_log" = "Nitrate",
                                "b_summerT" = "Temperature",
                                "b_NHD_RdDensWs" = "Road Density",
                                "b_log_width" = "log(Width)",
                                "b_Dam_binary1" = "Dam")) + #top
    theme_bw())

# Save out this figure.
# ggsave(y3_fig,
#        filename = "figures/teton_fall22/brms_yield_nuts_030123.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

# Plot conditional effects of covariates.

y3_n <- conditional_effects(y3, effects = "no3_log")

# Create new dataframe
y3n_df <- y3_n$no3_log %>%
  # and calculate true yield values
  # checked to be sure these default probs are 0.95
  mutate(yield = 10^`estimate__`,
         loweryield = 10^`lower__`,
         upperyield = 10^`upper__`)

(plot_y3n <- ggplot(y3n_df, aes(x = 10^no3_log, y = yield)) +
    geom_line(color = "black", linewidth = 1) +
    geom_ribbon(aes(ymin = loweryield, ymax = upperyield),
                alpha = 0.25) +
    geom_point(data = dat_yield_brms3, 
               aes(x = 10^no3_log, y = 10^log_yield),
               alpha = 0.2) +
    labs(x = expression(Mean~Nitrate~(mg/L~NO[3]-N)),
         y = expression(a[max])) +
    scale_x_log10() +
    ylim(0, 17) +
    theme_bw())

y3_p <- conditional_effects(y3, effects = "p_log")

# Create new dataframe
y3p_df <- y3_p$p_log %>%
  # and calculate true yield values
  # checked to be sure these default probs are 0.95
  mutate(yield = 10^`estimate__`,
         loweryield = 10^`lower__`,
         upperyield = 10^`upper__`)

(plot_y3p <- ggplot(y3p_df, aes(x = 10^p_log, y = yield)) +
    geom_line(color = "black", linewidth = 1) +
    geom_ribbon(aes(ymin = loweryield, ymax = upperyield),
                alpha = 0.25) +
    geom_point(data = dat_yield_brms3, 
               aes(x = 10^p_log, y = 10^log_yield),
               alpha = 0.2) +
    labs(x = expression(Mean~Dissolved~Phosphorus~(mg/L~P)),
         y = expression(a[max])) +
    scale_x_log10() +
    ylim(0, 17) +
    theme_bw())

y3_w <- conditional_effects(y3, effects = "log_width")

# Create new dataframe
y3w_df <- y3_w$log_width %>%
  # and calculate true Qc:Q2yr values
  mutate(yield = 10^`estimate__`,
         loweryield = 10^`lower__`,
         upperyield = 10^`upper__`)

(plot_y3w <- ggplot(y3w_df, aes(x = 10^log_width, y = yield)) +
    geom_line(color = "black", linewidth = 1) +
    geom_ribbon(aes(ymin = loweryield, ymax = upperyield),
                alpha = 0.25) +
    geom_point(data = dat_yield_brms3, 
               aes(x = 10^log_width, y = 10^log_yield),
               alpha = 0.2) +
    scale_x_log10()+
    labs(x = "River Width (m)",
         y = expression(a[max])) +
    ylim(0, 17) +
    theme_bw())

y3_t <- conditional_effects(y3, effects = "summerT")

# Create new dataframe
y3t_df <- y3_t$summerT %>%
  # and calculate true Qc:Q2yr values
  mutate(yield = 10^`estimate__`,
         loweryield = 10^`lower__`,
         upperyield = 10^`upper__`)

(plot_y3t <- ggplot(y3t_df, aes(x = summerT, y = yield)) +
    geom_line(color = "black", linewidth = 1) +
    geom_ribbon(aes(ymin = loweryield, ymax = upperyield),
                alpha = 0.25) +
    geom_point(data = dat_yield_brms3, 
               aes(x = summerT, y = 10^log_yield),
               alpha = 0.2) +
    #scale_x_log10()+
    labs(x = paste0("Mean Summer Temperature (", '\u00B0', "C)"),
         y = expression(a[max])) +
    ylim(0, 17) +
    theme_bw())

y3_rd <- conditional_effects(y3, effects = "NHD_RdDensWs")

# Create new dataframe
y3rd_df <- y3_rd$NHD_RdDensWs %>%
  # and calculate true Qc:Q2yr values
  mutate(yield = 10^`estimate__`,
         loweryield = 10^`lower__`,
         upperyield = 10^`upper__`)

(plot_y3rd <- ggplot(y3rd_df, aes(x = NHD_RdDensWs, y = yield)) +
    geom_line(color = "black", linewidth = 1) +
    geom_ribbon(aes(ymin = loweryield, ymax = upperyield),
                alpha = 0.25) +
    geom_point(data = dat_yield_brms3, 
               aes(x = NHD_RdDensWs, y = 10^log_yield),
               alpha = 0.2) +
    #scale_x_log10()+
    labs(x = expression(Road~Density~by~Watershed~(km/km^2)),
         y = expression(a[max])) +
    ylim(0, 17) +
    theme_bw())

y3_d <- conditional_effects(y3, effects = "Dam_binary")

# Create new dataframe
y3d_df <- y3_d$Dam_binary %>%
  # and calculate true yield values
  mutate(yield = 10^`estimate__`,
         loweryield = 10^`lower__`,
         upperyield = 10^`upper__`)

(plot_y3d <- ggplot(y3d_df, aes(x = Dam_binary, y = yield)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = loweryield, ymax = upperyield), width = 0.2) +
    geom_jitter(data = dat_yield_brms3, aes(x = Dam_binary, y = 10^log_yield),
                alpha = 0.1, width = 0.1) +
    labs(x = "Likelihood of Interference by Dams",
         y = expression(a[max])) +
    scale_x_discrete(labels = c("5-50%", "100%")) +
    ylim(0, 17) +
    theme_bw())

# Now, let's combine the above using patchwork.
(fig_cond_yield_nuts <- (y3_fig | (plot_y3d / plot_y3t) | (plot_y3n / plot_y3p) | (plot_y3rd / plot_y3w)) +
    plot_annotation(tag_levels = 'A'))

# And export.
# ggsave(fig_cond_yield_nuts,
#        filename = "figures/teton_fall22/brms_yield_cond_nuts_032323.jpg",
#        width = 40,
#        height = 20,
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
  dplyr::select(site_name,
         c_med, Qc_Q2yr, cvQ, meanGPP:NHD_AREASQKM, 
         NHD_RdDensCat:NHD_PctImp2011Ws, 
         Canal:width_med)

# Join with precip data.
dat_Qc_trim <- left_join(dat_Qc_trim, site_precip, 
                         by = c("site_name" = "SiteID"))

# Join with HUC2 data.
dat_Qc_trim <- left_join(dat_Qc_trim, site_HUC2)

# And visualize the relationships with Qc:Qy2r values.
QcQ2_covs <- ggpairs(dat_Qc_trim %>% 
                       dplyr::select(-site_name) %>%
                       dplyr::select(-NHD_RdDensCat,
                                     -NHD_PctImp2011Cat)) # trim for space

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
# to remove GPP.

# (3) Remaining Pearson's correlation values are below 0.5.

# Log transform and edit necessary covariates.
dat_Qc_trim <- dat_Qc_trim %>%
  mutate(GPP_log = log10(meanGPP),
         area_log = log10(NHD_AREASQKM),
         width_log = log10(width_med)) %>%
  # Also creating a new categorical dam column to model by.
  mutate(Dam_binary = factor(case_when(
    Dam %in% c("50", "80", "95") ~ "0", # Potential
    Dam == "0" ~ "1", # Certain
    TRUE ~ NA)))

# Notes on model structure:

# Qc:Q2yr ~ size + roads + dams + 1 | HUC2

# I am going to include the following covariates as representatives of the
# corresponding environmental factors:

# width - stream size
# road density - terrestrial development
# dam - aquatic development
# HUC2 - to account for random effect of geography

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
# precipitation - land and water connectivity metric but poor data availability

# One on one plots for covariates of interest vs. QcQ2yr.
hist(dat_Qc_trim$Qc_Q2yr)

# Going to log transform QcQ2yr too.
dat_Qc_trim <- dat_Qc_trim %>%
  mutate(logQcQ2 = log10(Qc_Q2yr))

plot(logQcQ2 ~ width_log, data = dat_Qc_trim)
plot(logQcQ2 ~ NHD_RdDensWs, data = dat_Qc_trim)
plot(logQcQ2 ~ Dam_binary, data = dat_Qc_trim)
plot(logQcQ2 ~ huc2_id, data = dat_Qc_trim)
hist(dat_Qc_trim$logQcQ2)

# Ok, and making the final dataset with which to build models.
dat_Qc_brms <- dat_Qc_trim %>%
  dplyr::select(logQcQ2, width_log, NHD_RdDensWs, Dam_binary, huc2_id)

##### Step 1: Create multi-level model.

# Remove NAs since STAN doesn't play nicely with them.
mq1 <- dat_Qc_brms %>%
  drop_na(NHD_RdDensWs, Dam_binary, width_log, huc2_id)

q1 <- brm(logQcQ2 ~ NHD_RdDensWs + Dam_binary + width_log + (1|huc2_id), 
          data = mq1, family = gaussian())
# assumes 4 chains and 2000 iterations (1000 warm-up)

# Export for safekeeping.
#saveRDS(q1, "data_posthoc_modelfits/qcq2_brms_030123.rds")

##### Step 2: Examine model outputs.

summary(q1)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept       -0.24      0.18    -0.60     0.11 1.00     2831     1834
# NHD_RdDensWs    -0.01      0.02    -0.05     0.02 1.00     3358     2540
# Dam_binary1     -0.06      0.09    -0.23     0.11 1.00     4436     3214
# width_log        0.02      0.11    -0.20     0.24 1.00     3063     2060

# Well, great convergence still, but nothing looks significant.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged well, so let's look at chain
# mixing and posterior distributions.
plot(q1)

# Chains all appear well-mixed, but let's also check things in shinystan.
launch_shinystan(q1)

# No divergent transitions :)

# Finally, examine to be sure no n_eff are < 0.1 or 10%
mcmc_plot(q1, type = "neff")

##### Step 4: Examine model relationships for each predictor.

plot(conditional_effects(q1, effects = "NHD_RdDensWs"))
plot(conditional_effects(q1, effects = "Dam_binary"))
plot(conditional_effects(q1, effects = "width_log"))

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
mq1$obs <- c(1:length(mq1$logQcQ2))

q1.1 <- brm(logQcQ2 ~ NHD_RdDensWs + Dam_binary + width_log + 
              (1|huc2_id) + (1|obs), 
            data = mq1, family = gaussian())
# 322 divergent transitions EEE!

# Compare with original model using leave-one-out approximation.
loo(q1, q1.1)

# Model comparisons:
#     elpd_diff se_diff
# q1.1   0.0       0.0  
# q1    -4.6       1.1 

# Higher expected log posterior density (elpd) values = better fit.
# So, in this case model accounting for overdispersion (q1.1) fits better.
# But there are 67 problematic observations and hundreds of divergent
# transitions so sticking with the original model (q1).

##### Step 6: Plot the results.

get_variables(q1)

# b_Intercept refers to global mean
# r_huc2_id[] are the offsets from that mean for each condition

(q_fig <- mcmc_plot(q1, variable = c("b_width_log", "b_NHD_RdDensWs", 
                                     "b_Dam_binary1"),
                    #type = "intervals",
                    point_est = "median", # default = "median"
                    prob = 0.66, # default = 0.5
                    prob_outer = 0.95) + # default = 0.9
    vline_at(v = 0) +
    labs(x = "Posterior Estimates",
         y = "Predictors") +
    scale_y_discrete(labels = c("b_NHD_RdDensWs" = "Road Density",
                                "b_Dam_binary1" = "Dam",
                                "b_width_log" = "log(Width)")) +
    theme_bw())

# Plot conditional effects of all covariates.
# Using code from here to make them ggplots:
# https://bookdown.org/content/4857/conditional-manatees.html#summary-bonus-conditional_effects

q_r <- conditional_effects(q1, effects = "NHD_RdDensWs")

# Create new dataframe
qr_df <- q_r$NHD_RdDensWs %>%
  # and calculate true Qc:Q2yr values
  mutate(QcQ2 = 10^`estimate__`,
         lowerQcQ2 = 10^`lower__`,
         upperQcQ2 = 10^`upper__`)

(plot_qr <- ggplot(qr_df, aes(x = NHD_RdDensWs, y = QcQ2)) +
    geom_line(color = "black", linewidth = 1) +
    geom_ribbon(aes(ymin = lowerQcQ2, ymax = upperQcQ2),
                alpha = 0.25) +
    geom_point(data = dat_Qc_trim, aes(x = NHD_RdDensWs, y = Qc_Q2yr),
                alpha = 0.2) +
    labs(x = expression(Road~Density~by~Watershed~(km/km^2)),
         y = expression(Q[c]:Q[2~yr])) +
    ylim(0, 1.2) +
    theme_bw())

q_d <- conditional_effects(q1, effects = "Dam_binary")

# Create new dataframe
qd_df <- q_d$Dam_binary %>%
  # and calculate true Qc:Q2yr values
  mutate(QcQ2 = 10^`estimate__`,
         lowerQcQ2 = 10^`lower__`,
         upperQcQ2 = 10^`upper__`)

(plot_qd <- ggplot(qd_df, aes(x = Dam_binary, y = QcQ2)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lowerQcQ2, ymax = upperQcQ2), 
                  width = 0.2) +
    geom_jitter(data = dat_Qc_trim, aes(x = Dam_binary, y = Qc_Q2yr),
                alpha = 0.1, width = 0.1) +
    labs(x = "Likelihood of Interference by Dams",
         y = expression(Q[c]:Q[2~yr])) +
    scale_x_discrete(labels = c("5-50%", "100%")) +
    ylim(0, 1.2) +
    theme_bw())

q_w <- conditional_effects(q1, effects = "width_log")

# Create new dataframe
qw_df <- q_w$width_log %>%
  # and calculate true Qc:Q2yr values
  mutate(QcQ2 = 10^`estimate__`,
         lowerQcQ2 = 10^`lower__`,
         upperQcQ2 = 10^`upper__`)

(plot_qw <- ggplot(qw_df, aes(x = 10^width_log, y = QcQ2)) +
    geom_line(color = "black", linewidth = 1) +
    geom_ribbon(aes(ymin = lowerQcQ2, ymax = upperQcQ2),
                alpha = 0.25) +
    geom_point(data = dat_Qc_trim, aes(x = 10^width_log, y = Qc_Q2yr),
               alpha = 0.2) +
    scale_x_log10()+
    labs(x = "River Width (m)",
         y = expression(Q[c]:Q[2~yr])) +
    ylim(0, 1.2) +
    theme_bw())

# Now, let's combine the above using patchwork.
(fig_cond_Qc <- (q_fig + plot_qd + plot_qr + plot_qw) +
    plot_layout(nrow = 1) +
    plot_annotation(tag_levels = 'A'))

# And export.
# ggsave(fig_cond_Qc,
#        filename = "figures/teton_fall22/brms_Qc_cond_032423.jpg",
#        width = 40,
#        height = 10,
#        units = "cm")

# Save out this figure.
# ggsave(q_fig,
#        filename = "figures/teton_fall22/brms_QcQ2_021623.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

# End of script.