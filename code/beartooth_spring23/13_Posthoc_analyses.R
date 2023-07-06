## Recovery of Stream Productivity following Disturbance
## Originally created: November 22, 2022
## Heili Lowman

#### READ ME ####

# The following set of scripts will walk through the steps necessary to
# prep and send data to Beartooth as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_XX" folders have been ignored using git.ignore.
# Links to the raw data sets are provided in the 01_NWIS_RiverSelection.R file.

# If you are accessing the code via GitHub, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running the 
# code below.

# This file will perform the post-hoc analyses to explore covariates that
# might explain the median parameter estimates from our model output.

# The following script was run using the Pinyon server at the University of 
# Nevada Reno for speed.

#### Setup ####

# Load necessary packages.
lapply(c("tidybayes", "brms", "tidyverse", "lubridate", 
         "data.table", "GGally",
         "multcomp", "patchwork", "bayesplot",
         "modelsummary", "here", "nlme","loo"), 
       require, character.only=T)

#### Data ####

# Import necessary datasets.

# First, the data for the maximum accrual (amax) models.
dat_amax <- readRDS("data_working/amax_covariates_152sites_070523.rds")

# Next, the data for the Qc:Q2yr models.
dat_Qc <- readRDS("data_working/Qc_covariates_138sites_070523.rds")

#### Model 1: Max. Algal Yield using 'brms' ####

# First, visualize the relationships with median yield values.
may_covs <- ggpairs(dat_amax %>% 
                       dplyr::select(yield_med, 
                                     cvQ:NHD_PctImp2011Ws,
                                     Canal:Phosphorus,
                                     huc_2:exc_y))

# ggsave(may_covs,
#        filename = "figures/beartooth_spring23/yield_covariates.jpg",
#        width = 50,
#        height = 50,
#        units = "cm")

# Some notes regarding these covariates.

# (1) I will transform variables to make parameter estimates comparable.

# (2) I will keep in mind correlations btw covariates found earlier.

# (3) The following jump out as potentially important (not including those
# data used to actually generate estimates): temperature, latitude, road
# density/pct impervious in watershed, maybe canals+dams??, width, P, and
# highest outliers appear in certain HUC2s, which may be an artefact of the
# geographic distribution of these sites - single digit HUC2s are on the 
# E coast for the most part.

# (4) Also checked to be sure there was a lack of correlation between # of days
# in a given record and the yield values.

# Proposed starting model structure:

# max algal yield/accrual ~ size + roads + dams + 
#                           temperature + exceedances/year + 1 | HUC2

# Will create a separate model adding in nutrients based on best model
# fit here (since records are far fewer).

hist(dat_amax$yield_med)

# Edit dataset.
dat_amax <- dat_amax %>%
  mutate(log_yield = log10(yield_med)) %>% # log-transform yield
  mutate(log_width = log10(width_med)) %>% # log-transform width
  # Also creating a new categorical dam column to model by.
  # same metric used in QC:Q2yr model below.
  mutate(Dam_binary = factor(case_when(
    Dam %in% c("50", "80", "95") ~ "0", # Potential = 5-50%
    Dam == "0" ~ "1", # Certain = 100%
    TRUE ~ NA)))

# One on one plots for covariates of interest vs. may.
plot(log_yield ~ meanTemp, data = dat_amax)
plot(log_yield ~ NHD_RdDensWs, data = dat_amax)
plot(log_yield ~ Dam_binary, data = dat_amax)
plot(log_yield ~ log_width, data = dat_amax)
plot(log_yield ~ huc_2, data = dat_amax)
plot(log_yield ~ exc_y, data = dat_amax)
hist(dat_amax$log_yield)

# Ok, and making the final dataset with which to build models
# where necessary variables have already been log-transformed and
# now just need to be scaled.
dat_amax_brms1 <- dat_amax %>%
  # assigning sites to be rownames so that we can re-identify and add HUC2
  # and Dams back in once we've scaled the remaining variables
  column_to_rownames(var = "site_name") %>%
  dplyr::select(log_yield, meanTemp, NHD_RdDensWs, log_width, exc_y)
  
dat_amax_brms1_scaled <- scale(dat_amax_brms1) # scale variables

# Pull sites back in so that we can match with HUC2 values.
dat_amax_brms <- rownames_to_column(as.data.frame(dat_amax_brms1_scaled), 
                                     var = "site_name")

dat_amax_Dam_HUC <- dat_amax %>%
  dplyr::select(site_name, Dam_binary, huc_2)

dat_amax_brms <- left_join(dat_amax_brms, dat_amax_Dam_HUC) %>%
  dplyr::select(log_yield, meanTemp, NHD_RdDensWs, 
                Dam_binary, log_width, exc_y, huc_2) %>%
  mutate(huc_2 = factor(huc_2))

##### Step 1: Create multi-level model.

a1 <- brm(log_yield ~ log_width + NHD_RdDensWs +
            Dam_binary + meanTemp + exc_y + (1|huc_2), 
         data = dat_amax_brms, family = gaussian())
# assumes 4 chains and 2000 iterations (1000 warm-up)
# Runs in about a minute on the server :)
# 8 sites omitted due to missing data

# Export for safekeeping.
# saveRDS(a1, "data_posthoc_modelfits/accrual_brms_070623.rds")

##### Step 2: Examine model outputs.

summary(a1)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept        0.23      0.14    -0.04     0.52 1.00     2896     2561
# log_width        0.52      0.10     0.33     0.72 1.00     4555     2928
# NHD_RdDensWs     0.03      0.10    -0.17     0.22 1.00     3880     3109
# Dam_binary1     -0.40      0.15    -0.70    -0.11 1.00     6863     3095
# meanTemp        -0.01      0.09    -0.18     0.17 1.00     4113     2911
# exc_y            0.16      0.07     0.02     0.31 1.00     5475     3101

# Well, for one, this is great convergence! All Rhat < 1.05.
# And at first glance, exceedance events, dam_binary1, and 
# stream width jump out as important.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged well, so let's look at chain
# mixing and posterior distributions.
plot(a1, variable = c("b_meanTemp", "b_NHD_RdDensWs",
                      "b_Dam_binary1", "b_log_width",
                      "b_exc_y"))

# Chains all appear well-mixed.

# Finally, examine to be sure no n_eff are < 0.1
mcmc_plot(a1, type = "neff")

##### Step 4: Examine model relationships for each predictor.

plot(conditional_effects(a1, effects = "log_width"))
plot(conditional_effects(a1, effects = "NHD_RdDensWs"))
plot(conditional_effects(a1, effects = "Dam_binary"))
plot(conditional_effects(a1, effects = "meanTemp"))
plot(conditional_effects(a1, effects = "exc_y"))

# Note, can investigate scenarios like effect1:effect2 here
# and it will automatically choose percentiles to predict.

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
dat_amax_brms$obs <- c(1:length(dat_amax_brms$log_yield))

a1.1 <- brm(log_yield ~ log_width + NHD_RdDensWs +
            Dam_binary + meanTemp + exc_y + (1|huc_2) + (1|obs), 
          data = dat_amax_brms, family = gaussian())
# 221 divergent transitions EEK!

# Compare with original model using leave-one-out approximation.
loo(a1, a1.1)

# Model comparisons:
#     elpd_diff se_diff
# a1.1   0.0       0.0  
# a1    -5.4       0.7 

# Higher expected log posterior density (elpd) values = better fit.
# So, in this case model accounting for overdispersion (a1.1) fits better.
# But there are 31 problematic observations and 200+ divergences,
# as well as high Rhat values, so my gut says the original model (a1) is better.

##### Step 6: Plot the results.

get_variables(a1)

# b_Intercept refers to global mean
# r_huc_2[] are the offsets from that mean for each condition

##### Figures #####

# Examine the data.
post_data <- mcmc_intervals_data(a1,
                         point_est = "median", # default = "median"
                         prob = 0.66, # default = 0.5
                         prob_outer = 0.95) # default = 0.9

View(post_data)

# Making custom plot to change color of each interval.
# Using core dataset "post_data" rather than canned function.
(a_fig_custom <- ggplot(post_data %>%
                        filter(parameter %in% c("b_log_width", "b_exc_y",
                                                "b_NHD_RdDensWs",
                                                "b_meanTemp", 
                                                "b_Dam_binary1")) %>%
                        mutate(par_f = factor(parameter, 
                                              levels = c("b_log_width",
                                                         "b_exc_y",
                                                         "b_NHD_RdDensWs",
                                                         "b_meanTemp",
                                                         "b_Dam_binary1"))), 
                        aes(x = m, y = par_f, color = par_f)) +
    geom_linerange(aes(xmin = ll, xmax = hh),
                  size = 2, alpha = 0.5) +
    geom_point(size = 5) +
    vline_at(v = 0) +
    scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
    labs(x = "Posterior Estimates",
         y = "Predictors") +
    scale_y_discrete(labels = c("b_log_width" = "Width",
                                "b_NHD_RdDensWs" = "Roads",
                                "b_Dam_binary1" = "Dam",
                                "b_meanTemp" = "Temperature",
                                "b_exc_y" = "Exceedances")) +
    theme_bw() +
    scale_color_manual(values = c("#4B8FF7", "#233D3F", "#233D3F", 
                                           "#4B8FF7", "#233D3F")) +
    theme(text = element_text(size = 10),
          legend.position = "none"))

# Can also use pars = c("^r_", "^b_", "^sd_") in place of variable phrasing
# to see all results.

# Default is type = "intervals".

# Also create conditional plots for each parameter.
# Plot conditional effects of all covariates.
# Using code from here to make them ggplots:
# https://bookdown.org/content/4857/conditional-manatees.html#summary-bonus-conditional_effects

# Also need to first "unscale" all values.

# Note the attributes of the originally scaled dataset.
center <- attr(dat_amax_brms1_scaled, "scaled:center")
scale <- attr(dat_amax_brms1_scaled, "scaled:scale")

# And use them to back-transform data below

###### Dams ######

# add_epred_draws() documentation here: https://mjskay.github.io/tidybayes/reference/add_predicted_draws.html, https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html

# Using some code from this resource to properly sample the binary:
# https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/#binary-effect

# Condition on dams alone.
a_d <- add_epred_draws(newdata = expand_grid(Dam_binary = c(0, 1),
                                    # hold remainder of covariates constant
                                    # choosing to predict for median/reference values
                                    NHD_RdDensWs = median(dat_amax_brms$NHD_RdDensWs,
                                                                   na.rm = TRUE),
                                    exc_y = median(dat_amax_brms$exc_y, 
                                                      na.rm = TRUE),
                                    log_width = median(dat_amax_brms$log_width, 
                                                       na.rm = TRUE),
                                    meanTemp = median(dat_amax_brms$meanTemp, 
                                                     na.rm = TRUE)),
                       object = a1,
                       re_formula = NA, # random effects not included due to 
                       # wanting to plot the global mean
                       ndraws = 100)

# Create new dataframe in the appropriate order as the scaled matrix.
ad_select <- a_d %>%
  ungroup() %>% # removing groups as imposed above
  dplyr::select(`.epred`, meanTemp, NHD_RdDensWs, log_width, exc_y) %>%
  rename("log_yield" = ".epred")

# And calculate true yield values
# Yield was scaled - Dam_binary was not
ad_descaled_data <- as.data.frame(t(t(ad_select) * scale + center)) %>%
  mutate(Dam_binary = factor(a_d$Dam_binary)) %>% # Add dams back in.
  mutate(draw = a_d$`.draw`) # Add draws back in.

# And plot all lines with original data points.
(plot_ad <- ggplot(ad_descaled_data, aes(x = Dam_binary, y = 10^log_yield)) +
    geom_point(size = 5, shape = 15, alpha = 0.2,  color = "#233D3F") +
    geom_jitter(data = dat_amax %>%
                  drop_na(Dam_binary), aes(x = Dam_binary, y = 10^log_yield),
                size = 3, alpha = 0.3, width = 0.1, color = "#233D3F") +
    labs(x = "Likelihood of Dam Interference",
         y = expression(a[max])) +
    scale_x_discrete(labels = c("5-50%", "100%")) +
    scale_y_log10() +
    theme_bw() +
    theme(text = element_text(size = 10)))

###### Temperature ######

# not plotting a spaghetti plot since there is no effect of temp on accrual 

(plot_at <- ggplot(dat_amax, aes(x = meanTemp, y = 10^log_yield)) +
    # Plot original unscaled data.
    geom_point(size = 3, alpha = 0.4, color = "#4B8FF7") +
    # And label things correctly.
    labs(x = paste0("Mean Summer Temperature (", '\u00B0', "C)"),
         y = expression(a[max])) +
    scale_y_log10() +
    theme_bw() +
    theme(text = element_text(size = 10)))

###### Roads ######

# not plotting a spaghetti plot since there is no effect of roads on accrual 

(plot_ard <- ggplot(dat_amax, aes(x = NHD_RdDensWs, y = 10^log_yield)) +
    geom_point(size = 3, alpha = 0.3, color = "#233D3F") +
    scale_y_log10()+
    labs(x = expression(Road~Density~(km/km^2)),
         y = expression(a[max])) +
    theme_bw() +
    theme(text = element_text(size = 10)))

###### Exceedances ######

# Condition on exceedances alone.
a_e <- add_epred_draws(newdata = expand_grid(exc_y = modelr::seq_range(dat_amax_brms$exc_y, n = 100),
                                  # hold remainder of covariates constant
                                  # choosing to predict for median/reference values
                                  NHD_RdDensWs = median(dat_amax_brms$NHD_RdDensWs,
                                                                   na.rm = TRUE),
                                  Dam_binary = c(0),
                                  log_width = median(dat_amax_brms$log_width, 
                                                     na.rm = TRUE),
                                  meanTemp = median(dat_amax_brms$meanTemp, 
                                                          na.rm = TRUE)),
                       object = a1,
                       re_formula = NA,
                       ndraws = 100)

# Create new dataframe in the appropriate order.
ae_select <- a_e %>%
  ungroup() %>% # removing groups as imposed above
  dplyr::select(`.epred`, meanTemp, NHD_RdDensWs, log_width, exc_y) %>%
  rename("log_yield" = ".epred")

# And calculate true yield values
ae_descaled_data <- as.data.frame(t(t(ae_select) * scale + center)) %>%
  mutate(draw = a_e$`.draw`) # Add draws back in.

# And plot all lines with original data points.
(plot_ae <- ggplot(ae_descaled_data, aes(x = exc_y, y = 10^log_yield)) +
    geom_line(aes(y = 10^log_yield, group = draw), alpha = 0.2, color = "#233D3F") +
    geom_point(data = dat_amax, aes(x = exc_y, y = 10^log_yield),
               size = 3, alpha = 0.3, color = "#233D3F") +
    scale_y_log10() +
    labs(x = expression(Mean~Annual~Exceedances~of~Q[c]),
         y = expression(a[max])) +
    theme_bw() +
    theme(text = element_text(size = 10)))

###### Width ######

# Condition on river width alone.
a_w <- add_epred_draws(newdata = expand_grid(log_width = modelr::seq_range(dat_amax_brms$log_width, n = 100),
                                  # hold remainder of covariates constant
                                  # choosing to predict for median/reference values
                                  NHD_RdDensWs = median(dat_amax_brms$NHD_RdDensWs,
                                                                   na.rm = TRUE),
                                  Dam_binary = c(0),
                                  meanTemp = median(dat_amax_brms$meanTemp, 
                                                          na.rm = TRUE),
                                  exc_y = median(dat_amax_brms$exc_y, na.rm = TRUE)),
                       object = a1,
                       re_formula = NA, # random effects not included due to global mean
                       ndraws = 100)

# Create new dataframe in the appropriate order to match scaled matrix.
aw_select <- a_w %>%
  ungroup() %>% # removing groups as imposed above
  dplyr::select(`.epred`, meanTemp, NHD_RdDensWs, log_width, exc_y) %>%
  rename("log_yield" = ".epred")

# And calculate true yield values
aw_descaled_data <- as.data.frame(t(t(aw_select) * scale + center)) %>%
  mutate(draw = a_w$`.draw`) # Add draws back in.

# And plot all lines with original data points.
(plot_aw <- ggplot(aw_descaled_data, aes(x = 10^log_width, y = 10^log_yield)) +
    geom_line(aes(y = 10^log_yield, group = draw), alpha = 0.2, color = "#4B8FF7") +
    geom_point(data = dat_amax, aes(x = 10^log_width, y = 10^log_yield),
               size = 3, alpha = 0.4, color = "#4B8FF7") +
    scale_y_log10() +
    scale_x_log10() +
    labs(x = "River Width (m)",
         y = expression(a[max])) +
    theme_bw() +
    theme(text = element_text(size = 10)))

###### Combined ######

# Now, to combine the above using patchwork.
(fig_cond_amax <- a_fig_custom + plot_ad + plot_at +
   plot_ard + plot_ae + plot_aw +
   plot_annotation(tag_levels = 'A'))

# And export.
# ggsave(fig_cond_amax,
#        filename = "figures/beartooth_spring23/brms_amax_cond_070623.jpg",
#        width = 22,
#        height = 13,
#        units = "cm")

##### Nutrients #####

# And build separate model for nutrients.
# Proposed starting model structure:

# yield ~ temp + roads + dams + width + exceedances + NO3 + P + 1 | HUC2

# Need to log transform nutrients.
dat_amax <- dat_amax %>%
  mutate(no3_log = log10(Nitrate),
         p_log = log10(Phosphorus))

# One on one plots for covariates of interest vs. may.
plot(log_yield ~ no3_log, data = dat_amax)
plot(log_yield ~ p_log, data = dat_amax)

# Ok, and making the final dataset with which to build models
dat_amax_brms2 <- dat_amax %>%
  # assigning sites to be rownames so that we can re-identify and add HUC2
  # back in once we've scaled the remaining variables
  column_to_rownames(var = "site_name") %>%
  dplyr::select(log_yield, meanTemp, NHD_RdDensWs,
                log_width, exc_y, no3_log, p_log)

dat_amax_brms2_scaled <- scale(dat_amax_brms2)

# Pull sites back in so that we can match with HUC2 values.
dat_amax_nuts_brms <- rownames_to_column(as.data.frame(dat_amax_brms2_scaled), 
                                     var = "site_name")

dat_amax_nuts_brms <- left_join(dat_amax_nuts_brms, dat_amax_Dam_HUC) %>%
  dplyr::select(log_yield, meanTemp, NHD_RdDensWs, 
                Dam_binary, log_width, exc_y,
                no3_log, p_log, huc_2) %>%
  mutate(huc_2 = factor(huc_2))

# And need to drop NAs.
dat_amax_nuts_brms <- dat_amax_nuts_brms %>%
  drop_na(no3_log, p_log) # 56 observations *sigh*

##### Step 1: Create multi-level model.

a2 <- brm(log_yield ~ meanTemp + NHD_RdDensWs + Dam_binary +
            log_width + exc_y + no3_log + p_log + (1|huc_2), 
          data = dat_amax_nuts_brms, family = gaussian(),
          iter = 2000) # 0 divergent transitions

# Export model fit for safekeeping.
# saveRDS(a2, "data_posthoc_modelfits/accrual_nuts_brms_070623.rds")

##### Step 2: Examine model outputs.

summary(a2)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept        0.05      0.17    -0.31     0.38 1.00     2463     2131
# meanTemp         0.07      0.14    -0.19     0.35 1.00     2847     2415
# NHD_RdDensWs    -0.15      0.14    -0.42     0.12 1.00     2539     2686
# Dam_binary1     -0.26      0.22    -0.70     0.19 1.00     2972     2564
# log_width        0.44      0.13     0.18     0.70 1.00     2605     2971
# exc_y            0.09      0.09    -0.09     0.27 1.00     3674     2708
# no3_log          0.07      0.16    -0.24     0.38 1.00     2922     2695
# p_log            0.16      0.14    -0.10     0.43 1.00     2432     2350

# Well, great convergence again! All Rhat < 1.05.
# And stream width again jumps out. But everything else seems to have
# been swamped out by the addition of nutrients.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged ok, so let's look at chain
# mixing and posterior distributions.
plot(a2, variable = c("b_meanTemp", "b_NHD_RdDensWs",
                      "b_Dam_binary1", "b_log_width",
                      "b_exc_y", "b_no3_log", "b_p_log"))

# Finally, examine to be sure no n_eff are < 0.1
mcmc_plot(a2, type = "neff")

##### Step 4: Examine model relationships for each predictor.

plot(conditional_effects(a2, effects = "meanTemp"))
plot(conditional_effects(a2, effects = "NHD_RdDensWs"))
plot(conditional_effects(a2, effects = "Dam_binary"))
plot(conditional_effects(a2, effects = "log_width"))
plot(conditional_effects(a2, effects = "exc_y"))
plot(conditional_effects(a2, effects = "no3_log"))
plot(conditional_effects(a2, effects = "p_log"))

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
dat_amax_nuts_brms$obs <- c(1:length(dat_amax_nuts_brms$log_yield))

a2.1 <- brm(log_yield ~ meanTemp + NHD_RdDensWs + Dam_binary +
              log_width + exc_y + no3_log + p_log +
              (1|huc_2) + (1|obs), 
            data = dat_amax_nuts_brms, family = gaussian())
# 154 divergent transitions eeeegads

# Compare with original model using leave-one-out approximation.
loo(a2, a2.1)

# Model comparisons:
#     elpd_diff se_diff
# a2.1   0.0       0.0  
# a2    -6.7       0.8

# Higher expected log posterior density (elpd) values = better fit.
# So, in this case model accounting for overdispersion (a2.1) fits better.
# But using the logic I employed above, I'm sticking with the original model 
# since it had FAR fewer divergent transitions.

##### Step 6: Plot the results.

get_variables(a2)

# b_Intercept refers to global mean
# r_huc2_id[] are the offsets from that mean for each condition

# Note the attributes of the originally scaled dataset.
center2 <- attr(dat_amax_brms2_scaled, "scaled:center")
scale2 <- attr(dat_amax_brms2_scaled, "scaled:scale")

###### Figures ######

# Examine the data.
post_data2 <- mcmc_intervals_data(a2,
                                 point_est = "median", # default = "median"
                                 prob = 0.66, # default = 0.5
                                 prob_outer = 0.95) # default = 0.9

View(post_data2)

# Making custom plot to change color of each interval.
# Using core dataset "post_data2" rather than canned function.
(a2_fig_custom <- ggplot(post_data2 %>%
                          filter(parameter %in% c("b_log_width", "b_exc_y",
                                                  "b_NHD_RdDensWs",
                                                  "b_meanTemp", "b_Dam_binary1",
                                                  "b_no3_log", "b_p_log")) %>%
                          mutate(par_f = factor(parameter, 
                                                levels = c("b_log_width",
                                                           "b_p_log",
                                                           "b_exc_y",
                                                           "b_meanTemp",
                                                           "b_no3_log",
                                                           "b_NHD_RdDensWs",
                                                           "b_Dam_binary1"))), 
                        aes(x = m, y = par_f, color = par_f)) +
    geom_linerange(aes(xmin = ll, xmax = hh),
                   size = 2, alpha = 0.5) +
    geom_point(size = 5) +
    vline_at(v = 0) +
    scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
    labs(x = "Posterior Estimates",
         y = "Predictors") +
    scale_y_discrete(labels = c("b_log_width" = "Width",
                                "b_NHD_RdDensWs" = "Roads",
                                "b_Dam_binary1" = "Dam",
                                "b_meanTemp" = "Temperature",
                                "b_exc_y" = "Exceedances",
                                "b_no3_log" = "Nitrate",
                                "b_p_log" = "Phosphorus")) +
    theme_bw() +
    scale_color_manual(values = c("#4B8FF7", "#4B8FF7", "#233D3F", "#4B8FF7",
                                  "#4B8FF7", "#233D3F", "#233D3F")) +
    theme(text = element_text(size = 10),
    legend.position = "none"))

###### Nitrate #######

# not plotting a spaghetti plot since there is no effect of NO3 on accrual 

(plot_a2n <- ggplot(dat_amax, aes(x = 10^no3_log, y = 10^log_yield)) +
    geom_point(size = 3, alpha = 0.4, color = "#4B8FF7") +
    labs(x = expression(Mean~Nitrate~(mg/L~NO[3]-N)),
         y = expression(a[max])) +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw() +
    theme(text = element_text(size = 10)))

####### Phosphorus #######

# not plotting a spaghetti plot since there is no effect of P on accrual 

(plot_a2p <- ggplot(dat_amax, aes(x = 10^p_log, y = 10^log_yield)) +
    geom_point(size = 3, alpha = 0.4, color = "#4B8FF7") +
    labs(x = expression(Mean~Phosphorus~(mg/L~P)),
         y = expression(a[max])) +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw() +
    theme(text = element_text(size = 10)))

####### Combined #######

# Now, let's combine the above using patchwork.
(fig_cond_amax_nuts <- a2_fig_custom + plot_a2n + plot_a2p +
    plot_annotation(tag_levels = 'A'))

# And export.
# ggsave(fig_cond_amax_nuts,
#        filename = "figures/beartooth_spring23/brms_amax_cond_nuts_070623.jpg",
#        width = 22,
#        height = 7,
#        units = "cm")

#### Model 2: Qc:Q2yr ####

# Log transform and edit necessary covariates.
dat_Qc <- dat_Qc %>%
  mutate(Qc_Q2yr = Qc/RI_2yr_Q_cms)

# And visualize the relationships with Qc:Qy2r values.
QcQ2_covs <- ggpairs(dat_Qc %>% 
                       dplyr::select(Qc_Q2yr, 
                                     cvQ:NHD_PctImp2011Ws,
                                     Canal:Phosphorus,
                                     huc_2:exc_y)) # trim for space

# ggsave(QcQ2_covs,
#        filename = "figures/beartooth_spring23/QcQ2_covariates.jpg",
#        width = 50,
#        height = 50,
#        units = "cm")

# Some notes regarding these covariates.

# (1) Road density and impervious cover appear tightly correlated (0.824), so
# I should probably only include one of these in the final model build again.

# (2) Stream width and meanGPP also appear correlated (0.618), more strongly
# than they did in the rmax covariate exploration. So, makes double sense
# to remove GPP.

# (3) The following jump out as potentially important (not including those
# data used to actually generate estimates): longitude ??

# Notes on model structure:

# Qc:Q2yr ~ width + roads + dams + 1 | HUC2

# One on one plots for covariates of interest vs. QcQ2yr.
hist(dat_Qc$Qc_Q2yr)

# Going to log transform QcQ2yr too.
dat_Qc <- dat_Qc %>%
  mutate(logQcQ2 = log10(Qc_Q2yr),
         log_width = log10(width_med)) %>%
  # Also creating the new categorical dam column to model by.
  mutate(Dam_binary = factor(case_when(
    Dam %in% c("50", "80", "95") ~ "0", # Potential
    Dam == "0" ~ "1", # Certain
    TRUE ~ NA)))

plot(logQcQ2 ~ log_width, data = dat_Qc)
plot(logQcQ2 ~ NHD_RdDensWs, data = dat_Qc)
plot(logQcQ2 ~ Dam_binary, data = dat_Qc)
plot(logQcQ2 ~ huc_2, data = dat_Qc)
hist(dat_Qc$logQcQ2)

# Ok, and making the final dataset with which to build models
# where necessary variables have already been log-transformed and
# now just need to be scaled.
dat_Qc_brms1 <- dat_Qc %>%
  # assigning sites to be rownames so that we can re-identify and add HUC2
  # back in once we've scaled the remaining variables
  column_to_rownames(var = "site_name") %>%
  dplyr::select(logQcQ2, NHD_RdDensWs, log_width)

dat_Qc_brms1_scaled <- scale(dat_Qc_brms1)

# Pull sites back in so that we can match with HUC2 values.
dat_Qc_brms <- rownames_to_column(as.data.frame(dat_Qc_brms1_scaled), 
                                  var = "site_name")

dat_Qc_Dam_HUC <- dat_Qc %>%
  dplyr::select(site_name, Dam_binary, huc_2)

dat_Qc_brms <- left_join(dat_Qc_brms, dat_Qc_Dam_HUC) %>%
  dplyr::select(logQcQ2, NHD_RdDensWs, 
                Dam_binary, log_width, huc_2) %>%
  mutate(huc_2 = factor(huc_2))

##### Step 1: Create multi-level model.

# Remove NAs since STAN doesn't play nicely with them.
dat_Qc_brms_noNA <- dat_Qc_brms %>%
  drop_na(logQcQ2)

q1 <- brm(logQcQ2 ~ NHD_RdDensWs + Dam_binary + log_width + (1|huc_2), 
          data = dat_Qc_brms_noNA, family = gaussian())
# assumes 4 chains and 2000 iterations (1000 warm-up)

# Export for safekeeping.
# saveRDS(q1, "data_posthoc_modelfits/qcq2_brms_051123.rds")

##### Step 2: Examine model outputs.

summary(q1)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept        0.20      0.21    -0.21     0.61 1.00     1772     1976
# NHD_RdDensWs    -0.06      0.11    -0.29     0.16 1.00     3110     2962
# Dam_binary1     -0.22      0.18    -0.57     0.13 1.00     5089     2785
# log_width        0.05      0.11    -0.18     0.27 1.00     3182     2896

# Well, great convergence still, but nothing looks significant.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged well, so let's look at chain
# mixing and posterior distributions.
plot(q1)

# Finally, examine to be sure no n_eff are < 0.1 or 10%
mcmc_plot(q1, type = "neff")

##### Step 4: Examine model relationships for each predictor.

plot(conditional_effects(q1, effects = "NHD_RdDensWs"))
plot(conditional_effects(q1, effects = "Dam_binary"))
plot(conditional_effects(q1, effects = "log_width"))

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
dat_Qc_brms_noNA$obs <- c(1:length(dat_Qc_brms_noNA$logQcQ2))

q1.1 <- brm(logQcQ2 ~ NHD_RdDensWs + Dam_binary + log_width + 
              (1|huc_2) + (1|obs), 
            data = dat_Qc_brms_noNA, family = gaussian())
# 150 divergent transitions EEE!

# Compare with original model using leave-one-out approximation.
loo(q1, q1.1)

# Model comparisons:
#     elpd_diff se_diff
# q1.1   0.0       0.0  
# q1   -22.0       0.9 

# Higher expected log posterior density (elpd) values = better fit.
# So, in this case model accounting for overdispersion (q1.1) fits better.
# But there are 44 problematic observations and 150 divergent
# transitions so sticking with the original model (q1).

##### Step 6: Plot the results.

get_variables(q1)

# b_Intercept refers to global mean
# r_huc2_id[] are the offsets from that mean for each condition

# Note the attributes of the originally scaled dataset.
center3 <- attr(dat_Qc_brms1_scaled, "scaled:center")
scale3 <- attr(dat_Qc_brms1_scaled, "scaled:scale")

###### Figures ######

# Examine the data.
post_data3 <- mcmc_intervals_data(q1,
                                  point_est = "median", # default = "median"
                                  prob = 0.66, # default = 0.5
                                  prob_outer = 0.95) # default = 0.9

View(post_data3)

# Making custom plot to change color of each interval.
# Using core dataset "post_data" rather than canned function.
(q1_fig_custom <- ggplot(post_data3 %>%
                           filter(parameter %in% c("b_log_width",
                                                   "b_NHD_RdDensWs",
                                                   "b_Dam_binary1")) %>%
                           mutate(par_f = factor(parameter, 
                                                 levels = c("b_log_width",
                                                            "b_NHD_RdDensWs",
                                                            "b_Dam_binary1"))), 
                         aes(x = m, y = par_f, color = par_f)) +
    geom_linerange(aes(xmin = ll, xmax = hh),
                   size = 2, alpha = 0.5) +
    geom_point(size = 5) +
    vline_at(v = 0) +
    scale_x_continuous(breaks = c(-0.4,-0.2, 0, 0.2)) +
    labs(x = "Posterior Estimates",
         y = "Predictors") +
    scale_y_discrete(labels = c("b_log_width" = "Width",
                                "b_NHD_RdDensWs" = "Roads",
                                "b_Dam_binary1" = "Dam")) +
    theme_bw() +
    scale_color_manual(values = c("#4B8FF7", "#233D3F", "#233D3F")) +
    theme(text = element_text(size = 10),
    legend.position = "none"))

####### Roads #######

# not plotting a spaghetti plot since there is no effect of roads on Qc 

(plot_qr <- ggplot(dat_Qc, aes(x = NHD_RdDensWs, y = 10^logQcQ2)) +
    geom_point(size = 3, alpha = 0.3, color = "#233D3F") +
    labs(x = expression(Road~Density~(km/km^2)),
         y = expression(Q[c]:Q[2~yr])) +
    scale_y_log10() +
    theme_bw() +
    theme(text = element_text(size = 10)))

####### Dams #######

# not plotting overlapping draws since there is no effect of dams on Qc

(plot_qd <- ggplot(dat_Qc %>%
                     drop_na(Dam_binary), aes(x = Dam_binary, y = 10^logQcQ2)) +
    geom_jitter(size = 3, alpha = 0.3, width = 0.1, color = "#233D3F") +
    labs(x = "Likelihood of Dam Interference",
         y = expression(Q[c]:Q[2~yr])) +
    scale_x_discrete(labels = c("5-50%", "100%")) +
    scale_y_log10() +
    theme_bw() +
    theme(text = element_text(size = 10)))

####### Size #######

# not plotting a spaghetti plot since there is no effect of width on Qc 

(plot_qw <- ggplot(dat_Qc, aes(x = 10^log_width, y = 10^logQcQ2)) +
    geom_point(size = 3, alpha = 0.4, color = "#4B8FF7") +
    scale_y_log10() +
    scale_x_log10()+
    labs(x = "River Width (m)",
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw() +
    theme(text = element_text(size = 10)))

####### Combined #######

# Now, let's combine the above using patchwork.
(fig_cond_Qc <- (q1_fig_custom + plot_qd + plot_qr + plot_qw) +
    plot_layout(nrow = 2) +
    plot_annotation(tag_levels = 'A'))

# And export.
# ggsave(fig_cond_Qc,
#        filename = "figures/beartooth_spring23/brms_Qc_cond_051223.jpg",
#        width = 14,
#        height = 13,
#        units = "cm")

# End of script.
