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

# (3) Then, summer temp and latitude appear tightly correlated (0.666).
# They do provide somewhat different information, so I will include them as
# an interaction term. (Thinking about latitude as a temp gradient while
# longitude might be a kind of aridity gradient.)

# (4) After that, NO3 and PO4 appear the most tightly correlated (0.548).
# Those aren't duplicates, so I will include them as an interaction term.

# (5) Remaining Pearson's correlation values are below 0.5.

# Log transform necessary covariates.

dat_rmax_trim <- dat_rmax_trim %>%
  mutate(r_log = log10(r_med),
         GPP_log = log10(meanGPP),
         area_log = log10(NHD_AREASQKM),
         width_log = log10(width_med),
         no3_log = log10(Nitrate),
         po4_log = log10(Orthophosphate))

# Notes on model structure:

# I am going to include the following covariates as representatives of the
# corresponding environmental factors:

# cvQ - flow
# GPP - biological productivity (function of flow & light)*
# summer light - light at the stream surface during greatest canopy
# summer temp - water temperature (decoupled from air temperature)
# longitude - a rough location/aridity index?
# order - stream size**
# watershed area - stream size**
# width - stream size**
# road density - terrestrial development
# dam - aquatic development

# * Models will explore with and without GPP since GPP ~ f(flow, light).
# ** Models will explore each of the size indicators separately.

# The following covariates have been removed for the following reasons:

# latitude - too closely correlated with summer temperature
# % impervious land cover - too closely correlated with road density
# canal - another metric of terr/aq development but felt duplicative
# NO3/PO4 - 50% of sites have no data, so these will be examined separately

##### Stream Size #####

# Build initial set of models to investigate size covariates.

str(dat_rmax_trim) # keeping order, canal, and dam data as categorical

lm1_rmax <- lm(r_med ~ cvQ + GPP_log + summerL + summerT +
                 Lon_WGS84 + Order + NHD_RdDensCat + Dam, 
               data = dat_rmax_trim)

lm2_rmax <- lm(r_med ~ cvQ + GPP_log + summerL + summerT +
                 Lon_WGS84 + area_log + NHD_RdDensCat + Dam, 
               data = dat_rmax_trim)

lm3_rmax <- lm(r_med ~ cvQ + GPP_log + summerL + summerT +
                 Lon_WGS84 + width_log + NHD_RdDensCat + Dam, 
               data = dat_rmax_trim)

# Examine the outputs.

summary(lm1_rmax)
summary(lm2_rmax)
summary(lm3_rmax)

# Examine the coefficients.

lm1_rmax_tidy <- broom::tidy(lm1_rmax) # using order
View(lm1_rmax_tidy) # GPP, Lon, Temp, cvQ, dam (80%), 6th order (p<0.05)

lm2_rmax_tidy <- broom::tidy(lm2_rmax) # using area
View(lm2_rmax_tidy) # GPP, Lon, Temp, dam (80%, 95%), cvQ (p<0.05)

lm3_rmax_tidy <- broom::tidy(lm3_rmax) # using width
View(lm3_rmax_tidy) # GPP, Lon, Temp, dam(80%), cvQ (p<0.05)

# Examine model fit.

lm1_rmax_fit <- broom::glance(lm1_rmax) # order
View(lm1_rmax_fit) # adj R2 = 0.42, sigma = 0.08, p < 0.0001, nobs = 152

lm2_rmax_fit <- broom::glance(lm2_rmax) # area
View(lm2_rmax_fit) # adj R2 = 0.41, sigma = 0.08, p < 0.0001, nobs = 152

lm3_rmax_fit <- broom::glance(lm3_rmax) # width
View(lm3_rmax_fit) #adj R2 = 0.41, sigma = 0.08, p < 0.0001, nobs = 152

# Examine model diagnostics.
plot(lm1_rmax) # site 30 does appear to be an outlier here
plot(lm2_rmax) # looks a bit better
plot(lm3_rmax) # also fine

# The first residuals plots all look a bit curvy, so may think about log-transforming r values.

# Compare models.
aic1 <- AIC(lm1_rmax) # -317.45
aic2 <- AIC(lm2_rmax) # -320.39
aic3 <- AIC(lm3_rmax) # -316.98

# Moving forward with width as the indicator of stream size.

##### GPP #####

# Build second set of models to investigate influence of GPP (with/without).

# lm3_rmax <- lm(r_med ~ cvQ + GPP_log + summerL + summerT +
#                  Lon_WGS84 + width_log + NHD_RdDensCat + Dam, 
#                data = dat_rmax_trim)

lm4_rmax <- lm(r_med ~ cvQ + summerL + summerT +
                 Lon_WGS84 + width_log + NHD_RdDensCat + Dam, 
               data = dat_rmax_trim)

# Examine the outputs.

summary(lm3_rmax)
summary(lm4_rmax) # size jumps in importance.

# Examine the coefficients.

# with GPP, using width
View(lm3_rmax_tidy) # GPP, Lon, Temp, dam(80%), cvQ (p<0.05)

lm4_rmax_tidy <- broom::tidy(lm4_rmax) # without GPP, using width
View(lm4_rmax_tidy) # Width, Temp, Lon (p<0.05)

# Examine model fit.

# with GPP, width
View(lm3_rmax_fit) #adj R2 = 0.41, sigma = 0.08, p < 0.0001, nobs = 152

lm4_rmax_fit <- broom::glance(lm4_rmax) # without GPP, width
View(lm4_rmax_fit) #adj R2 = 0.15, sigma = 0.10, p = 0.0002, nobs = 151

# Examine model diagnostics.
plot(lm3_rmax) # fine
plot(lm4_rmax) # also fine - less curvy on first residuals pane

# Compare models.
aic3 # -316.98
aic4 <- AIC(lm4_rmax) # -262.86 EEK!

# Moving forward with GPP included in the model.

# Building one final model with rmax values on the log scale.

# lm3_rmax <- lm(r_med ~ cvQ + GPP_log + summerL + summerT +
#                  Lon_WGS84 + width_log + NHD_RdDensCat + Dam, 
#                data = dat_rmax_trim)

lm5_rmax <- lm(r_log ~ cvQ + GPP_log + summerL + summerT +
                 Lon_WGS84 + width_log + NHD_RdDensCat + Dam, 
               data = dat_rmax_trim)

# Examine the outputs.

summary(lm3_rmax)
summary(lm5_rmax) # more covariates emerge as important

# Examine the coefficients.

# with GPP, using width
View(lm3_rmax_tidy) # GPP, Lon, Temp, dam(80%), cvQ (p<0.05)

lm5_rmax_tidy <- broom::tidy(lm5_rmax) # with GPP, using width, log(rmax)
View(lm5_rmax_tidy) # GPP, Lon, cvQ, dam (80%), Rd Dens (p<0.05)

# Examine model fit.

# with GPP, width, rmax
View(lm3_rmax_fit) #adj R2 = 0.41, sigma = 0.08, p < 0.0001, nobs = 152

lm5_rmax_fit <- broom::glance(lm5_rmax) # with GPP, using width, log(rmax)
View(lm5_rmax_fit) #adj R2 = 0.51, sigma = 0.24, p < 0.0001, nobs = 151

# Examine model diagnostics.
plot(lm3_rmax) # curvy first pane that I'm trying to address
plot(lm5_rmax) # better residuals pane and slightly better QQ plot

# Compare models.
aic3 # -316.98
aic5 <- AIC(lm5_rmax) # 13.87 Oh dear!

# Moving forward with lm3_rmax.

# Export table of current results.
write_csv(lm3_rmax_tidy, "data_working/lm_rmax.csv")

##### Nutrients #####

# And build separate model for nutrients.
lmnuts_rmax <- lm(r_med ~ no3_log*po4_log, 
                  data = dat_rmax_trim)

# Examine the outputs.
summary(lmnuts_rmax)

# Examine the coefficients.

lmnuts_rmax_tidy <- broom::tidy(lmnuts_rmax)
View(lmnuts_rmax_tidy) # neither NO3 nor PO4 is significant

# Examine model fit.

lmnuts_rmax_fit <- broom::glance(lmnuts_rmax)
View(lmnuts_rmax_fit) # adj R2 = 0.09, sigma = 0.09, p = 0.01, nobs = 89

# Examine model diagnostics.
plot(lmnuts_rmax) # site 135 does *nearly* appear to be an outlier

#### Model 2: rmax vs. cvQ residuals ####

# Fit linear model between rmax median values and cvQ.
fit1 <- lm(r_med ~ cvQ, data = dat_rmax_trim)

dat_rmax_trim$residuals <- residuals(fit1)

# Building final model with residuals instead of rmax values.

lm1_resids <- lm(residuals ~ GPP_log + summerL + summerT +
                 Lon_WGS84 + width_log + NHD_RdDensCat + Dam,
               data = dat_rmax_trim)

# Examine the outputs.

summary(lm1_resids) # similar covariates emerge as important

# Examine the coefficients.

lm1_resids_tidy <- broom::tidy(lm1_resids) # with GPP, without cvQ, using residuals
View(lm1_resids_tidy) # GPP, Lon, dam (80%) (p<0.05)

# Examine model fit.

lm1_resids_fit <- broom::glance(lm1_resids) # with GPP, without cvQ, using residuals
View(lm1_resids_fit) #adj R2 = 0.34, sigma = 0.08, p < 0.0001, nobs = 151

# Examine model diagnostics.
plot(lm1_resids) # first pane still a bit curvy, same outliers as other lms (30, 91, 133)

# Compare models.
aic1r <- AIC(lm1_resids) # -304.46 Oh dear!

# Not sure this adds information since the significant parameters stay the same.

#### Model 3: Qc:Q2yr ####

# End of script.