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

# and the raw data used in model fitting.
dat_in <- readRDS("data_working/list_182sites_Qmaxnorm_allSL.rds")

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

# For ease, making the original list into a df.
dat_in_df <- map_df(dat_in, ~as.data.frame(.x), .id="site_name")

# And count years in each record at each site.
dat_years <- dat_in_df %>%
  group_by(site_name) %>%
  summarize(years = n_distinct(year)) %>%
  ungroup()

#### Model 1: Max. Algal Yield using 'brms' ####

# First, need to bind yield estimates with remaining data
# Using calculation based on eq. 7b from Scheuerell 2016 = `yield_med2`
# Combining with datasets from above.
dat_yield_combo1 <- left_join(dat_yield, site_precip, 
                             by = c("site_name" = "SiteID"))

dat_yield_combo2 <- left_join(dat_yield_combo1, site_HUC2) # HUC2

dat_yield_combo3 <- left_join(dat_yield_combo2, site_P) # dissolved P

dat_yield_combo4 <- left_join(dat_yield_combo3, dat_exc) # exceedances

dat_yield_combo <- left_join(dat_yield_combo4, dat_years) # # of years

# Add column to standardize exceedance events by number of years in
# a timeseries at a given site.
dat_yield_combo <- dat_yield_combo %>%
  mutate(exc_ev_y = total_exc_events/years)

# And visualize the relationships with median yield values.
may_covs <- ggpairs(dat_yield_combo %>% 
                       dplyr::select(yield_med2, cvQ:Orthophosphate,
                                     Phosphorus,
                                     pre_mm_cyr:huc2_id,
                                     exc_ev_y))

# ggsave(may_covs,
#        filename = "figures/teton_fall22/yield_covariates.jpg",
#        width = 50,
#        height = 50,
#        units = "cm")

# Some notes regarding these covariates.

# (1) I will transform variables to make parameter estimates comparable.

# (2) I will keep in mind correlations btw covariates found earlier.

# (3) The following jump out as potentially important (not including those
# data used to actually generate estimates): temperature, latitude, road
# density/pct impervious in watershed, empirical coefficient `a` of
# width estimation equation, maybe canals+dams??, width, precip, and highest
# outliers appear in certain HUC2s.

# Proposed starting model structure:

# max algal yield/accrual ~ size + roads + dams + temperature + exceedances
# + 1 | HUC2

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
plot(log_yield ~ exc_ev_y, data = dat_yield_combo)
hist(dat_yield_combo$log_yield)

# Ok, and making the final dataset with which to build models
# where necessary variables have already been log-transformed and
# now just need to be scaled.
dat_yield_brms1 <- dat_yield_combo %>%
  # assigning sites to be rownames so that we can re-identify and add HUC2
  # back in once we've scaled the remaining variables
  column_to_rownames(var = "site_name") %>%
  dplyr::select(log_yield, summerT, NHD_RdDensWs, log_width, exc_ev_y)
  
dat_yield_brms1 <- scale(dat_yield_brms1)

# Pull sites back in so that we can match with HUC2 values.
dat_yield_brms <- rownames_to_column(as.data.frame(dat_yield_brms1), 
                                     var = "site_name")

dat_sites_HUCs <- dat_yield_combo %>%
  dplyr::select(site_name, Dam_binary, huc2_id)

dat_yield_brms <- left_join(dat_yield_brms, dat_sites_HUCs) %>%
  dplyr::select(log_yield, summerT, NHD_RdDensWs, 
                Dam_binary, log_width, exc_ev_y, huc2_id) %>%
  mutate(huc2_id = factor(huc2_id))

##### Step 1: Create multi-level model.

y1 <- brm(log_yield ~ log_width + NHD_RdDensWs +
            Dam_binary + summerT + exc_ev_y + (1|huc2_id), 
         data = dat_yield_brms, family = gaussian())
# assumes 4 chains and 2000 iterations (1000 warm-up)
# started at 2:13pm - finished at 2:14pm on the server :)

# Export for safekeeping.
# saveRDS(y1, "data_posthoc_modelfits/accrual_brms_032923.rds")

##### Step 2: Examine model outputs.

summary(y1)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept        0.26      0.15    -0.03     0.57 1.00     2467     2172
# log_width        0.51      0.10     0.33     0.70 1.00     3359     3124
# NHD_RdDensWs     0.10      0.10    -0.10     0.28 1.00     3111     2537
# Dam_binary1     -0.40      0.15    -0.71    -0.10 1.00     4784     2984
# summerT         -0.25      0.08    -0.40    -0.08 1.00     3965     2565
# exc_ev_y         0.24      0.07     0.11     0.39 1.00     5804     2903

# Well, for one, this is great convergence! All Rhat < 1.05.
# And at first glance, exceedance events, temperature, dam_binary1, and 
# stream width jump out as important.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged well, so let's look at chain
# mixing and posterior distributions.
plot(y1, variable = c("b_summerT", "b_NHD_RdDensWs",
                      "b_Dam_binary1", "b_log_width",
                      "b_exc_ev_y"))

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
plot(conditional_effects(y1, effects = "exc_ev_y"))

# Note, can investigate scenarios like effect1:effect2 here
# and it will automatically choose percentiles to predict.

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
dat_yield_brms$obs <- c(1:length(dat_yield_brms$log_yield))

y1.1 <- brm(log_yield ~ log_width + NHD_RdDensWs +
            Dam_binary + summerT + exc_ev_y + (1|huc2_id) + (1|obs), 
          data = dat_yield_brms, family = gaussian())
# 353 divergent transitions EEK!

# Compare with original model using leave-one-out approximation.
loo(y1, y1.1)

# Model comparisons:
#     elpd_diff se_diff
# y1.1   0.0       0.0  
# y1    -6.1       1.0 

# Higher expected log posterior density (elpd) values = better fit.
# So, in this case model accounting for overdispersion (y1.1) fits better.
# But there are 30 problematic observations...
# And a few hundred divergences and poor Rhat values, so my gut says the
# original model (y1) has the better fit.

##### Step 6: Plot the results.

get_variables(y1)

# b_Intercept refers to global mean
# r_huc2_id[] are the offsets from that mean for each condition

##### Figures #####

# Examine the data being pulled by the function below
post_data <- mcmc_intervals_data(y1,
                         point_est = "median", # default = "median"
                         prob = 0.66, # default = 0.5
                         prob_outer = 0.95) # default = 0.9

View(post_data)

(y_fig <- mcmc_plot(y1, variable = c("b_log_width", "b_exc_ev_y",
                                     "b_NHD_RdDensWs",
                                     "b_summerT", "b_Dam_binary1"),
      #type = "intervals",
      point_est = "median", # default = "median"
      prob = 0.66, # default = 0.5
      prob_outer = 0.95) + # default = 0.9
    vline_at(v = 0) +
    labs(x = "Posterior Estimates",
         y = "Predictors") +
    scale_y_discrete(labels = c("b_log_width" = "Width",
                                "b_NHD_RdDensWs" = "Roads",
                                "b_Dam_binary1" = "Dam",
                                "b_summerT" = "Temperature",
                                "b_exc_ev_y" = "Exceedances")) +
    theme_bw())

# Save out this figure.
# ggsave(y_fig,
#        filename = "figures/teton_fall22/brms_yield_033123.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

# Making custom plot to change color of each interval.
# Using core dataset "post_data" rather than canned function.
(y_fig_custom <- ggplot(post_data %>%
                        filter(parameter %in% c("b_log_width", "b_exc_ev_y",
                                                "b_NHD_RdDensWs",
                                                "b_summerT", "b_Dam_binary1")) %>%
                        mutate(par_f = factor(parameter, 
                                              levels = c("b_log_width",
                                                         "b_exc_ev_y",
                                                         "b_NHD_RdDensWs",
                                                         "b_summerT",
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
                                "b_summerT" = "Temperature",
                                "b_exc_ev_y" = "Exceedances")) +
    theme_bw() +
    scale_color_manual(values = c("#4B8FF7", "#233D3F", "#233D3F", 
                                           "#4B8FF7", "#233D3F")) +
    theme(text = element_text(size = 24),
          legend.position = "none"))

# Save out this figure.
# ggsave(y_fig_custom,
#        filename = "figures/teton_fall22/brms_yield_custom_041723.jpg",
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

# Also need to first "unscale" all values.

# Note the attributes of the originally scaled dataset.
center <- attr(dat_yield_brms1, "scaled:center")
scale <- attr(dat_yield_brms1, "scaled:scale")

# And use them to back-transform data below

###### Dams ######

y_d <- conditional_effects(y1, effects = "Dam_binary")

# Create new dataframe
yd_df <- y_d$Dam_binary

# Yield was scaled - Dam_binary was not.
yd_select <- yd_df %>%
  dplyr::select(`estimate__`, summerT, NHD_RdDensWs, log_width, exc_ev_y) %>%
  rename("log_yield" = "estimate__")

# And calculate true yield values
yd_descaled_data <- as.data.frame(t(t(yd_select) * scale + center)) %>%
  mutate(Dam_binary = yd_df$Dam_binary)

# Also, need to do this for each of the 95% CIs, but the order of the
# de-scaling matters, so doing this twice more with each of the
# intervals as the first column.

# 2.5% lower interval
yd_select25 <- yd_df %>%
  dplyr::select(`lower__`, summerT, NHD_RdDensWs, log_width, exc_ev_y) %>%
  rename("lower_yield" = "lower__")

yd_descaled_data25 <- as.data.frame(t(t(yd_select25) * scale + center)) %>%
  mutate(Dam_binary = yd_df$Dam_binary) %>%
  dplyr::select(Dam_binary, lower_yield)

# 97.5% lower interval
yd_select975 <- yd_df %>%
  dplyr::select(`upper__`, summerT, NHD_RdDensWs, log_width, exc_ev_y) %>%
  rename("upper_yield" = "upper__")

yd_descaled_data975 <- as.data.frame(t(t(yd_select975) * scale + center)) %>%
  mutate(Dam_binary = yd_df$Dam_binary) %>%
  dplyr::select(Dam_binary, upper_yield)

yd_descaled_data <- left_join(yd_descaled_data, yd_descaled_data25)
yd_descaled_data <- left_join(yd_descaled_data, yd_descaled_data975)

(plot_yd <- ggplot(yd_descaled_data, aes(x = Dam_binary, y = 10^log_yield)) +
    geom_point(size = 5, color = "#233D3F") +
    geom_linerange(aes(ymin = 10^lower_yield, ymax = 10^upper_yield), 
                  size = 2, alpha = 0.7,
                  color = "#233D3F") +
    geom_jitter(data = dat_yield_combo, aes(x = Dam_binary, y = yield_med2),
                alpha = 0.3, width = 0.1, color = "#233D3F") +
    labs(x = "Likelihood of Interference by Dams",
         y = expression(a[max])) +
    scale_x_discrete(labels = c("5-50%", "100%")) +
    scale_y_log10() +
    theme_bw() +
    theme(text = element_text(size = 24)))

###### Temperature ######

# Working but still scaled example of spaghetti plot with 100 posterior means
# add_epred_draws() documentation here: https://mjskay.github.io/tidybayes/reference/add_predicted_draws.html, https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html
(testing <- add_epred_draws(newdata = expand_grid(summerT = modelr::seq_range(dat_yield_brms$summerT, 
                                                                              n = 100),
                                                  # hold remainder of covariates constant
                                                  # choosing to predict for median/reference values
                                                  NHD_RdDensWs = median(dat_yield_brms$NHD_RdDensWs,
                                                                        na.rm = TRUE),
                                                  Dam_binary = c(0),
                                                  log_width = median(dat_yield_brms$log_width,
                                                                     na.rm = TRUE),
                                                  exc_ev_y = median(dat_yield_brms$exc_ev_y,
                                                                    na.rm = TRUE)),
                            object = y1,
                            re_formula = NA, # random effects not included due to global mean
                            ndraws = 100) %>% # draw 100 different lines
   ggplot(aes(x = summerT, y = log_yield)) +
   geom_line(aes(y = .epred, group = paste(.draw)), alpha = 0.1, color = "#4B8FF7") +
   geom_point(data = dat_yield_brms, size = 3, alpha = 0.5, color = "#4B8FF7") +
   theme_bw())

# Save out this figure.
# ggsave(testing,
#        filename = "figures/teton_fall22/brms_temp_predlines_041723.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

# Condition on temperature alone.
y_t <- add_epred_draws(newdata = expand_grid(summerT = modelr::seq_range(dat_yield_brms$summerT, n = 100),
                       # hold remainder of covariates constant
                       # choosing to predict for median/reference values
                       NHD_RdDensWs = median(dat_yield_brms$NHD_RdDensWs,
                                                                   na.rm = TRUE),
                       Dam_binary = c(0),
                       log_width = median(dat_yield_brms$log_width, na.rm = TRUE),
                       exc_ev_y = median(dat_yield_brms$exc_ev_y, na.rm = TRUE)),
                       object = y1,
                       re_formula = NA, # random effects not included due to global mean
                       ndraws = 100)

# Create new dataframe in the appropriate order.
yt_select <- y_t %>%
  ungroup() %>% # removing groups as imposed above
  dplyr::select(`.epred`, summerT, NHD_RdDensWs, log_width, exc_ev_y) %>%
  rename("log_yield" = ".epred")

# And calculate true yield values
yt_descaled_data <- as.data.frame(t(t(yt_select) * scale + center)) %>%
  mutate(draw = y_t$`.draw`) # Add draws back in.

# And plot all lines with original data points.
(plot_yt <- ggplot(yt_descaled_data, aes(x = summerT, y = 10^log_yield)) +
    # Plot 100 lines.
    geom_line(aes(y = 10^log_yield, group = draw), alpha = 0.2, color = "#4B8FF7") +
    # Plot original unscaled data.
    geom_point(data = dat_yield_combo, aes(x = summerT, y = 10^log_yield),
               alpha = 0.4, color = "#4B8FF7") +
    # And label things correctly.
    labs(x = paste0("Mean Summer Temperature (", '\u00B0', "C)"),
         y = expression(a[max])) +
    scale_y_log10() +
    theme_bw() +
    theme(text = element_text(size = 24)))

###### Roads ######

y_rd <- conditional_effects(y1, effects = "NHD_RdDensWs")

# Create new dataframe
yrd_df <- y_rd$NHD_RdDensWs

yrd_select <- yrd_df %>%
  dplyr::select(`estimate__`, summerT, `effect1__`, log_width, exc_ev_y) %>%
  rename("log_yield" = "estimate__",
         "NHD_RdDensWs" = "effect1__")

# And calculate true yield values
yrd_descaled_data <- as.data.frame(t(t(yrd_select) * scale + center))

# Also, need to do this for each of the 95% CIs, but the order of the
# de-scaling matters, so doing this twice more with each of the
# intervals as the first column.

# 2.5% lower interval
yrd_select25 <- yrd_df %>%
  dplyr::select(`lower__`, summerT, `effect1__`, log_width, exc_ev_y) %>%
  rename("lower_yield" = "lower__",
         "NHD_RdDensWs" = "effect1__")

yrd_descaled_data25 <- as.data.frame(t(t(yrd_select25) * scale + center)) %>%
  dplyr::select(NHD_RdDensWs, lower_yield)

# 97.5% lower interval
yrd_select975 <- yrd_df %>%
  dplyr::select(`upper__`, summerT, `effect1__`, log_width, exc_ev_y) %>%
  rename("upper_yield" = "upper__",
         "NHD_RdDensWs" = "effect1__")

yrd_descaled_data975 <- as.data.frame(t(t(yrd_select975) * scale + center)) %>%
  dplyr::select(NHD_RdDensWs, upper_yield)

yrd_descaled_data <- left_join(yrd_descaled_data, yrd_descaled_data25)
yrd_descaled_data <- left_join(yrd_descaled_data, yrd_descaled_data975)

(plot_yrd <- ggplot(yrd_descaled_data, aes(x = NHD_RdDensWs, y = 10^log_yield)) +
    #geom_line(color = "black", linewidth = 1) +
    #geom_ribbon(aes(ymin = 10^lower_yield, ymax = 10^upper_yield),
    #            alpha = 0.25) +
    geom_point(data = dat_yield_combo, 
               aes(x = NHD_RdDensWs, y = yield_med2),
               alpha = 0.3, color = "#233D3F") +
    scale_y_log10()+
    labs(x = expression(Watershed~Road~Density~(km/km^2)),
         y = expression(a[max])) +
    theme_bw() +
    theme(text = element_text(size = 24)))

###### Exceedances ######

y_e <- conditional_effects(y1, effects = "exc_ev_y")

# Create new dataframe
ye_df <- y_e$exc_ev_y

ye_select <- ye_df %>%
  dplyr::select(`estimate__`, summerT, NHD_RdDensWs, log_width, `effect1__`) %>%
  rename("log_yield" = "estimate__",
         "exc_ev_y" = "effect1__")

# And calculate true yield values
ye_descaled_data <- as.data.frame(t(t(ye_select) * scale + center))

# Also, need to do this for each of the 95% CIs, but the order of the
# de-scaling matters, so doing this twice more with each of the
# intervals as the first column.

# 2.5% lower interval
ye_select25 <- ye_df %>%
  dplyr::select(`lower__`, summerT, NHD_RdDensWs, log_width, `effect1__`) %>%
  rename("lower_yield" = "lower__",
         "exc_ev_y" = "effect1__")

ye_descaled_data25 <- as.data.frame(t(t(ye_select25) * scale + center)) %>%
  dplyr::select(exc_ev_y, lower_yield)

# 97.5% lower interval
ye_select975 <- ye_df %>%
  dplyr::select(`upper__`, summerT, NHD_RdDensWs, log_width, `effect1__`) %>%
  rename("upper_yield" = "upper__",
         "exc_ev_y" = "effect1__")

ye_descaled_data975 <- as.data.frame(t(t(ye_select975) * scale + center)) %>%
  dplyr::select(exc_ev_y, upper_yield)

ye_descaled_data <- left_join(ye_descaled_data, ye_descaled_data25)
ye_descaled_data <- left_join(ye_descaled_data, ye_descaled_data975)

(plot_ye <- ggplot(ye_descaled_data, aes(x = exc_ev_y, y = 10^log_yield)) +
    geom_line(color = "#233D3F", linewidth = 1) +
    geom_ribbon(aes(ymin = 10^lower_yield, ymax = 10^upper_yield),
                alpha = 0.25, color = "#233D3F", fill = "#233D3F") +
    geom_point(data = dat_yield_combo, 
               aes(x = exc_ev_y, y = yield_med2),
               alpha = 0.3, color = "#233D3F") +
    scale_y_log10() +
    labs(x = expression(Annual~Exceedances~of~Q[c]),
         y = expression(a[max])) +
    theme_bw() +
    theme(text = element_text(size = 24)))

###### Size ######

y_w <- conditional_effects(y1, effects = "log_width")

# Create new dataframe
yw_df <- y_w$log_width

yw_select <- yw_df %>%
  dplyr::select(`estimate__`, summerT, NHD_RdDensWs, `effect1__`, exc_ev_y) %>%
  rename("log_yield" = "estimate__",
         "log_width" = "effect1__")

# And calculate true yield values
yw_descaled_data <- as.data.frame(t(t(yw_select) * scale + center))

# Also, need to do this for each of the 95% CIs, but the order of the
# de-scaling matters, so doing this twice more with each of the
# intervals as the first column.

# 2.5% lower interval
yw_select25 <- yw_df %>%
  dplyr::select(`lower__`, summerT, NHD_RdDensWs, `effect1__`, exc_ev_y) %>%
  rename("lower_yield" = "lower__",
         "log_width" = "effect1__")

yw_descaled_data25 <- as.data.frame(t(t(yw_select25) * scale + center)) %>%
  dplyr::select(log_width, lower_yield)

# 97.5% lower interval
yw_select975 <- yw_df %>%
  dplyr::select(`upper__`, summerT, NHD_RdDensWs, `effect1__`, exc_ev_y) %>%
  rename("upper_yield" = "upper__",
         "log_width" = "effect1__")

yw_descaled_data975 <- as.data.frame(t(t(yw_select975) * scale + center)) %>%
  dplyr::select(log_width, upper_yield)

yw_descaled_data <- left_join(yw_descaled_data, yw_descaled_data25)
yw_descaled_data <- left_join(yw_descaled_data, yw_descaled_data975)

(plot_yw <- ggplot(yw_descaled_data, aes(x = 10^log_width, y = 10^log_yield)) +
    geom_line(color = "#4B8FF7", linewidth = 1) +
    geom_ribbon(aes(ymin = 10^lower_yield, ymax = 10^upper_yield),
                alpha = 0.25, fill = "#4B8FF7", color = "#4B8FF7") +
    geom_point(data = dat_yield_combo, aes(x = width_med, y = yield_med2),
               alpha = 0.3, color = "#4B8FF7") +
    scale_y_log10() +
    scale_x_log10() +
    labs(x = "River Width (m)",
         y = expression(a[max])) +
    theme_bw() +
    theme(text = element_text(size = 24)))

###### Combined ######

# Now, let's combine the above using patchwork.
(fig_cond_yield <- y_fig_custom + plot_yd + plot_yt +
   plot_yrd + plot_ye + plot_yw +
   plot_annotation(tag_levels = 'A'))

# And export.
# ggsave(fig_cond_yield,
#        filename = "figures/teton_fall22/brms_yield_cond_041723.jpg",
#        width = 55,
#        height = 30,
#        units = "cm")

##### Nutrients #####

# And build separate model for nutrients.
# Proposed starting model structure:

# yield ~ temp + roads + dams + width + exceedances + NO3 + P + 1 | HUC2

# Need to log transform nutrients.
dat_yield_combo <- dat_yield_combo %>%
  mutate(no3_log = log10(Nitrate),
         p_log = log10(Phosphorus))

# One on one plots for covariates of interest vs. may.
plot(log_yield ~ no3_log, data = dat_yield_combo)
plot(log_yield ~ p_log, data = dat_yield_combo)

# Ok, and making the final dataset with which to build models
dat_yield_brms3 <- dat_yield_combo %>%
  # assigning sites to be rownames so that we can re-identify and add HUC2
  # back in once we've scaled the remaining variables
  column_to_rownames(var = "site_name") %>%
  dplyr::select(log_yield, summerT, NHD_RdDensWs,
                log_width, exc_ev_y, no3_log, p_log)

dat_yield_brms3 <- scale(dat_yield_brms3)

# Pull sites back in so that we can match with HUC2 values.
dat_yield_brms2 <- rownames_to_column(as.data.frame(dat_yield_brms3), 
                                     var = "site_name")

dat_sites_HUCs <- dat_yield_combo %>%
  dplyr::select(site_name, Dam_binary, huc2_id)

dat_yield_brms2 <- left_join(dat_yield_brms2, dat_sites_HUCs) %>%
  dplyr::select(log_yield, summerT, NHD_RdDensWs, 
                Dam_binary, log_width, exc_ev_y,
                no3_log, p_log, huc2_id) %>%
  mutate(huc2_id = factor(huc2_id))

# And need to drop NAs.
dat_yield_brms2 <- dat_yield_brms2 %>%
  drop_na(no3_log, p_log) # 60 observations *sigh*

##### Step 1: Create multi-level model.

y2 <- brm(log_yield ~ summerT + NHD_RdDensWs + Dam_binary +
            log_width + exc_ev_y + no3_log + p_log + (1|huc2_id), 
          data = dat_yield_brms2, family = gaussian(),
          iter = 2000) # 1 divergent transition

# Export model fit for safekeeping.
# saveRDS(y2, "data_posthoc_modelfits/accrual_nuts_brms_033023.rds")

##### Step 2: Examine model outputs.

summary(y2)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept        0.14      0.18    -0.21     0.48 1.00     1744     1462
# summerT          0.01      0.17    -0.31     0.34 1.00     2187     1772
# NHD_RdDensWs     0.01      0.15    -0.27     0.32 1.00     1204      908
# Dam_binary1     -0.36      0.21    -0.77     0.06 1.00     2548     1914
# log_width        0.50      0.14     0.22     0.78 1.00     2025     1476
# exc_ev_y         0.18      0.09     0.00     0.35 1.00     3698     2588
# no3_log         -0.06      0.16    -0.39     0.26 1.00     2011     2060
# p_log            0.18      0.14    -0.10     0.46 1.00     2079     2465

# Well, great convergence again! All Rhat < 1.05.
# And stream width again jumps out. But everything else seems to have
# been swamped out by the addition of nutrients.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged ok, so let's look at chain
# mixing and posterior distributions.
plot(y2, variable = c("b_summerT", "b_NHD_RdDensWs",
                      "b_Dam_binary1", "b_log_width",
                      "b_exc_ev_y",
                      "b_no3_log", "b_p_log"))

# Chains all appear well-mixed, but let's also check things in shinystan.
launch_shinystan(y2)
# 2 divergent transitions.

# Finally, examine to be sure no n_eff are < 0.1
mcmc_plot(y2, type = "neff")

##### Step 4: Examine model relationships for each predictor.

plot(conditional_effects(y2, effects = "summerT"))
plot(conditional_effects(y2, effects = "NHD_RdDensWs"))
plot(conditional_effects(y2, effects = "Dam_binary"))
plot(conditional_effects(y2, effects = "log_width"))
plot(conditional_effects(y2, effects = "exc_ev_y"))
plot(conditional_effects(y2, effects = "no3_log"))
plot(conditional_effects(y2, effects = "p_log"))

# Note, can investigate scenarios like effect1:effect2 here
# and it will automatically choose percentiles to predict.

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
dat_yield_brms2$obs <- c(1:length(dat_yield_brms2$log_yield))

y2.1 <- brm(log_yield ~ summerT + NHD_RdDensWs + Dam_binary +
              log_width + exc_ev_y + no3_log + p_log +
              (1|huc2_id) + (1|obs), 
            data = dat_yield_brms2, family = gaussian())
# 87 divergent transitions eeeegads

# Compare with original model using leave-one-out approximation.
loo(y2, y2.1)

# Model comparisons:
#     elpd_diff se_diff
# y2.1   0.0       0.0  
# y2   -12.0       1.2 

# Higher expected log posterior density (elpd) values = better fit.

# So, in this case model accounting for overdispersion (y2.1) fits better.
# But there are 19 problematic observations, and using the logic I
# employed above, I'm sticking with the original model since it
# had FAR fewer divergent transitions and better convergence (i.e., Rhat).

##### Step 6: Plot the results.

get_variables(y2)

# b_Intercept refers to global mean
# r_huc2_id[] are the offsets from that mean for each condition

# Note the attributes of the originally scaled dataset.
center2 <- attr(dat_yield_brms3, "scaled:center")
scale2 <- attr(dat_yield_brms3, "scaled:scale")

###### Figures ######

color_scheme_set("teal")

(y2_fig <- mcmc_plot(y2, variable = c("b_log_width", "b_exc_ev_y", "b_p_log",
                                      "b_summerT", "b_NHD_RdDensWs",
                                      "b_no3_log", "b_Dam_binary1"),
                     point_est = "median", # default = "median"
                     prob = 0.66, # default = 0.5
                     prob_outer = 0.95) + # default = 0.9
    vline_at(v = 0) +
    labs(x = "Posterior Estimates",
         y = "Predictors") +
    scale_y_discrete(labels = c("b_p_log" = "Phosphorus",
                                "b_no3_log" = "Nitrate",
                                "b_exc_ev_y" = "Exceedances",
                                "b_summerT" = "Temperature",
                                "b_NHD_RdDensWs" = "Roads",
                                "b_log_width" = "Width",
                                "b_Dam_binary1" = "Dam")) + #top
    theme_bw())

# Save out this figure.
# ggsave(y2_fig,
#        filename = "figures/teton_fall22/brms_yield_nuts_033123.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

####### Nitrate #######

# Plot conditional effects of additional nutrient covariates.

y2_n <- conditional_effects(y2, effects = "no3_log")

# Create new dataframe
y2n_df <- y2_n$no3_log 

y2n_select <- y2n_df %>%
  dplyr::select(`estimate__`, summerT, NHD_RdDensWs, log_width, 
                exc_ev_y, `effect1__`, p_log) %>%
  rename("log_yield" = "estimate__",
         "no3_log" = "effect1__")

# And calculate true yield values
y2n_descaled_data <- as.data.frame(t(t(y2n_select) * scale2 + center2))

# Also, need to do this for each of the 95% CIs, but the order of the
# de-scaling matters, so doing this twice more with each of the
# intervals as the first column.

# 2.5% lower interval
y2n_select25 <- y2n_df %>%
  dplyr::select(`lower__`, summerT, NHD_RdDensWs, log_width, 
                exc_ev_y, `effect1__`, p_log) %>%
  rename("lower_yield" = "lower__",
         "no3_log" = "effect1__")

y2n_descaled_data25 <- as.data.frame(t(t(y2n_select25) * scale2 + center2)) %>%
  dplyr::select(no3_log, lower_yield)

# 97.5% lower interval
y2n_select975 <- y2n_df %>%
  dplyr::select(`upper__`, summerT, NHD_RdDensWs, log_width, 
                exc_ev_y, `effect1__`, p_log) %>%
  rename("upper_yield" = "upper__",
         "no3_log" = "effect1__")

y2n_descaled_data975 <- as.data.frame(t(t(y2n_select975) * scale2 + center2)) %>%
  dplyr::select(no3_log, upper_yield)

y2n_descaled_data <- left_join(y2n_descaled_data, y2n_descaled_data25)
y2n_descaled_data <- left_join(y2n_descaled_data, y2n_descaled_data975)

(plot_y2n <- ggplot(y2n_descaled_data, aes(x = 10^no3_log, y = 10^log_yield)) +
    # geom_line(color = "#D46F10", linewidth = 1) +
    # geom_ribbon(aes(ymin = 10^lower_yield, ymax = 10^upper_yield),
    #             alpha = 0.25, color = "#D46F10", fill = "#D46F10") +
    geom_point(data = dat_yield_combo, 
               aes(x = 10^no3_log, y = 10^log_yield),
               alpha = 0.3, color = "#D46F10") +
    labs(x = expression(Mean~Nitrate~(mg/L~NO[3]-N)),
         y = expression(a[max])) +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw())

####### Phosphorus #######

y2_p <- conditional_effects(y2, effects = "p_log")

# Create new dataframe
y2p_df <- y2_p$p_log

y2p_select <- y2p_df %>%
  dplyr::select(`estimate__`, summerT, NHD_RdDensWs, log_width, 
                exc_ev_y, no3_log, `effect1__`) %>%
  rename("log_yield" = "estimate__",
         "p_log" = "effect1__")

# And calculate true yield values
y2p_descaled_data <- as.data.frame(t(t(y2p_select) * scale2 + center2))

# Also, need to do this for each of the 95% CIs, but the order of the
# de-scaling matters, so doing this twice more with each of the
# intervals as the first column.

# 2.5% lower interval
y2p_select25 <- y2p_df %>%
  dplyr::select(`lower__`, summerT, NHD_RdDensWs, log_width, 
                exc_ev_y, no3_log, `effect1__`) %>%
  rename("lower_yield" = "lower__",
         "p_log" = "effect1__")

y2p_descaled_data25 <- as.data.frame(t(t(y2p_select25) * scale2 + center2)) %>%
  dplyr::select(p_log, lower_yield)

# 97.5% lower interval
y2p_select975 <- y2p_df %>%
  dplyr::select(`upper__`, summerT, NHD_RdDensWs, log_width, 
                exc_ev_y, no3_log, `effect1__`) %>%
  rename("upper_yield" = "upper__",
         "p_log" = "effect1__")

y2p_descaled_data975 <- as.data.frame(t(t(y2p_select975) * scale2 + center2)) %>%
  dplyr::select(p_log, upper_yield)

y2p_descaled_data <- left_join(y2p_descaled_data, y2p_descaled_data25)
y2p_descaled_data <- left_join(y2p_descaled_data, y2p_descaled_data975)

(plot_y2p <- ggplot(y2p_descaled_data, aes(x = 10^p_log, y = 10^log_yield)) +
    # geom_line(color = "#D46F10", linewidth = 1) +
    # geom_ribbon(aes(ymin = 10^lower_yield, ymax = 10^upper_yield),
    #             alpha = 0.25, color = "#D46F10", fill = "#D46F10") +
    geom_point(data = dat_yield_combo, 
               aes(x = 10^p_log, y = 10^log_yield),
               alpha = 0.3, color = "#D46F10") +
    labs(x = expression(Mean~Dissolved~Phosphorus~(mg/L~P)),
         y = expression(a[max])) +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw())

####### Combined #######

# Now, let's combine the above using patchwork.
(fig_cond_yield_nuts <- y2_fig + plot_y2n + plot_y2p +
    plot_annotation(tag_levels = 'A'))

# And export.
# ggsave(fig_cond_yield_nuts,
#        filename = "figures/teton_fall22/brms_yield_cond_nuts_033123.jpg",
#        width = 40,
#        height = 10,
#        units = "cm")

# Also experimenting with a way to plot histograms/density plots
# over top of these intervals.

my_posterior2 <- as.matrix(y2)

(y2_fig2 <- mcmc_areas(my_posterior2,
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
# ggsave(y2_fig2,
#        filename = "figures/teton_fall22/brms_yield3_sig_021623.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

##### Without Outliers #####

# Following lab meeting, I'll do a quick check to see if removing the two
# highest outliers of accrual changes the results of the initial model build.

dat_yield_brms1.1 <- dat_yield_combo %>%
  # assigning sites to be rownames so that we can re-identify and add HUC2
  # back in once we've scaled the remaining variables
  column_to_rownames(var = "site_name") %>%
  dplyr::select(log_yield, summerT, NHD_RdDensWs, log_width, exc_ev_y) %>%
  filter(log_yield < 1.4)
  
dat_yield_brms1.1 <- scale(dat_yield_brms1.1)

# Pull sites back in so that we can match with HUC2 values.
dat_yield_brms0.1 <- rownames_to_column(as.data.frame(dat_yield_brms1.1), 
                                     var = "site_name")

dat_sites_HUCs <- dat_yield_combo %>%
  dplyr::select(site_name, Dam_binary, huc2_id)

dat_yield_brms0.1 <- left_join(dat_yield_brms0.1, dat_sites_HUCs) %>%
  dplyr::select(log_yield, summerT, NHD_RdDensWs, 
                Dam_binary, log_width, exc_ev_y, huc2_id) %>%
  mutate(huc2_id = factor(huc2_id))

# Fit model.

y1.1 <- brm(log_yield ~ log_width + NHD_RdDensWs +
            Dam_binary + summerT + exc_ev_y + (1|huc2_id), 
          data = dat_yield_brms0.1, family = gaussian())

# Export for safekeeping.
saveRDS(y1.1, "data_posthoc_modelfits/accrual_brms_no_outliers_033123.rds")

# Examine model outputs.
summary(y1.1)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept        0.32      0.16     0.02     0.65 1.00     2132     2386
# log_width        0.50      0.09     0.31     0.68 1.00     3654     3263
# NHD_RdDensWs     0.11      0.10    -0.09     0.30 1.00     3943     3245
# Dam_binary1     -0.49      0.15    -0.78    -0.20 1.00     5443     2971
# summerT         -0.25      0.09    -0.42    -0.09 1.00     3729     2756
# exc_ev_y         0.28      0.07     0.14     0.42 1.00     5654     3084

# And compare results figure.

color_scheme_set("teal")

(y0.1_fig <- mcmc_plot(y1.1, variable = c("b_log_width", "b_exc_ev_y",
                                     "b_NHD_RdDensWs",
                                     "b_summerT", "b_Dam_binary1"),
                    point_est = "median", # default = "median"
                    prob = 0.66, # default = 0.5
                    prob_outer = 0.95) + # default = 0.9
    vline_at(v = 0) +
    labs(x = "Posterior Estimates",
         y = "Predictors") +
    scale_y_discrete(labels = c("b_log_width" = "Width",
                                "b_NHD_RdDensWs" = "Roads",
                                "b_Dam_binary1" = "Dam",
                                "b_summerT" = "Temperature",
                                "b_exc_ev_y" = "Exceedances")) +
    theme_bw())

# Save out this figure.
# ggsave(y0.1_fig,
#        filename = "figures/teton_fall22/brms_yield_no_outliers_033123.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

#### Model 2: Qc:Q2yr ####

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

# Log transform and edit necessary covariates.
dat_Qc_trim <- dat_Qc_trim %>%
  mutate(GPP_log = log10(meanGPP),
         area_log = log10(NHD_AREASQKM),
         log_width = log10(width_med)) %>%
  # Also creating a new categorical dam column to model by.
  mutate(Dam_binary = factor(case_when(
    Dam %in% c("50", "80", "95") ~ "0", # Potential
    Dam == "0" ~ "1", # Certain
    TRUE ~ NA)))

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

plot(logQcQ2 ~ log_width, data = dat_Qc_trim)
plot(logQcQ2 ~ NHD_RdDensWs, data = dat_Qc_trim)
plot(logQcQ2 ~ Dam_binary, data = dat_Qc_trim)
plot(logQcQ2 ~ huc2_id, data = dat_Qc_trim)
hist(dat_Qc_trim$logQcQ2)

# Ok, and making the final dataset with which to build models
# where necessary variables have already been log-transformed and
# now just need to be scaled.
dat_Qc_brms1 <- dat_Qc_trim %>%
  # assigning sites to be rownames so that we can re-identify and add HUC2
  # back in once we've scaled the remaining variables
  column_to_rownames(var = "site_name") %>%
  dplyr::select(logQcQ2, NHD_RdDensWs, log_width)

dat_Qc_brms1 <- scale(dat_Qc_brms1)

# Pull sites back in so that we can match with HUC2 values.
dat_Qc_brms <- rownames_to_column(as.data.frame(dat_Qc_brms1), 
                                  var = "site_name")

dat_sites_HUCs2 <- dat_Qc_trim %>%
  dplyr::select(site_name, Dam_binary, huc2_id)

dat_Qc_brms <- left_join(dat_Qc_brms, dat_sites_HUCs2) %>%
  dplyr::select(logQcQ2, NHD_RdDensWs, 
                Dam_binary, log_width, huc2_id) %>%
  mutate(huc2_id = factor(huc2_id))

##### Step 1: Create multi-level model.

# Remove NAs since STAN doesn't play nicely with them.
mq1 <- dat_Qc_brms %>%
  drop_na()

q1 <- brm(logQcQ2 ~ NHD_RdDensWs + Dam_binary + log_width + (1|huc2_id), 
          data = mq1, family = gaussian())
# assumes 4 chains and 2000 iterations (1000 warm-up)

# Export for safekeeping.
#saveRDS(q1, "data_posthoc_modelfits/qcq2_brms_033123.rds")

##### Step 2: Examine model outputs.

summary(q1)

#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept        0.12      0.15    -0.17     0.42 1.00     2412     2265
# NHD_RdDensWs    -0.09      0.10    -0.30     0.12 1.00     3519     3184
# Dam_binary1     -0.12      0.19    -0.48     0.26 1.00     4211     2583
# log_width        0.02      0.11    -0.20     0.23 1.00     3225     2945

# Well, great convergence still, but nothing looks significant.

##### Step 3: Examine model diagnostics.

# Everything appears to have converged well, so let's look at chain
# mixing and posterior distributions.
plot(q1)

# Chains all appear well-mixed, but let's also check things in shinystan.
launch_shinystan(q1)

# Only 2 divergent transitions :)

# Finally, examine to be sure no n_eff are < 0.1 or 10%
mcmc_plot(q1, type = "neff")

##### Step 4: Examine model relationships for each predictor.

plot(conditional_effects(q1, effects = "NHD_RdDensWs"))
plot(conditional_effects(q1, effects = "Dam_binary"))
plot(conditional_effects(q1, effects = "log_width"))

##### Step 5: Investigate possible overdispersion.

# Add column denoting number of observations.
mq1$obs <- c(1:length(mq1$logQcQ2))

q1.1 <- brm(logQcQ2 ~ NHD_RdDensWs + Dam_binary + log_width + 
              (1|huc2_id) + (1|obs), 
            data = mq1, family = gaussian())
# 60 divergent transitions EEE!

# Compare with original model using leave-one-out approximation.
loo(q1, q1.1)

# Model comparisons:
#     elpd_diff se_diff
# q1.1   0.0       0.0  
# q1    -7.0       1.0 

# Higher expected log posterior density (elpd) values = better fit.
# So, in this case model accounting for overdispersion (q1.1) fits better.
# But there are 43 problematic observations and hundreds of divergent
# transitions so sticking with the original model (q1).

##### Step 6: Plot the results.

get_variables(q1)

# b_Intercept refers to global mean
# r_huc2_id[] are the offsets from that mean for each condition

# Note the attributes of the originally scaled dataset.
center3 <- attr(dat_Qc_brms1, "scaled:center")
scale3 <- attr(dat_Qc_brms1, "scaled:scale")

###### Figures ######

color_scheme_set("teal")
(q_fig <- mcmc_plot(q1, variable = c("b_log_width", "b_NHD_RdDensWs", 
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
                                "b_log_width" = "Width")) +
    theme_bw())

# Plot conditional effects of all covariates.
# Using code from here to make them ggplots:
# https://bookdown.org/content/4857/conditional-manatees.html#summary-bonus-conditional_effects

####### Roads #######

q_r <- conditional_effects(q1, effects = "NHD_RdDensWs")

# Create new dataframe
qr_df <- q_r$NHD_RdDensWs

qr_select <- qr_df %>%
  dplyr::select(`estimate__`, `effect1__`, log_width) %>%
  rename("logQcQ2" = "estimate__",
         "NHD_RdDensWs" = "effect1__")

# And calculate true yield values
qr_descaled_data <- as.data.frame(t(t(qr_select) * scale3 + center3))

# Also, need to do this for each of the 95% CIs, but the order of the
# de-scaling matters, so doing this twice more with each of the
# intervals as the first column.

# 2.5% lower interval
qr_select25 <- qr_df %>%
  dplyr::select(`lower__`, `effect1__`, log_width) %>%
  rename("lower_QcQ2" = "lower__",
         "NHD_RdDensWs" = "effect1__")

qr_descaled_data25 <- as.data.frame(t(t(qr_select25) * scale3 + center3)) %>%
  dplyr::select(NHD_RdDensWs, lower_QcQ2)

# 97.5% lower interval
qr_select975 <- qr_df %>%
  dplyr::select(`upper__`, `effect1__`, log_width) %>%
  rename("upper_QcQ2" = "upper__",
         "NHD_RdDensWs" = "effect1__")

qr_descaled_data975 <- as.data.frame(t(t(qr_select975) * scale3 + center3)) %>%
  dplyr::select(NHD_RdDensWs, upper_QcQ2)

qr_descaled_data <- left_join(qr_descaled_data, qr_descaled_data25)
qr_descaled_data <- left_join(qr_descaled_data, qr_descaled_data975)

(plot_qr <- ggplot(qr_descaled_data, aes(x = NHD_RdDensWs, y = 10^logQcQ2)) +
    # geom_line(color = "black", linewidth = 1) +
    # geom_ribbon(aes(ymin = lowerQcQ2, ymax = upperQcQ2),
    #             alpha = 0.25) +
    geom_point(data = dat_Qc_trim, aes(x = NHD_RdDensWs, y = Qc_Q2yr),
                alpha = 0.4, color = "#5A7ECB") +
    labs(x = expression(Road~Density~by~Watershed~(km/km^2)),
         y = expression(Q[c]:Q[2~yr])) +
    scale_y_log10() +
    theme_bw())

####### Dams #######

q_d <- conditional_effects(q1, effects = "Dam_binary")

# Create new dataframe
qd_df <- q_d$Dam_binary

# Qc:Q2 was scaled - Dam_binary was not.
qd_select <- qd_df %>%
  dplyr::select(`estimate__`, NHD_RdDensWs, log_width) %>%
  rename("logQcQ2" = "estimate__")

# And calculate true yield values
qd_descaled_data <- as.data.frame(t(t(qd_select) * scale3 + center3)) %>%
  mutate(Dam_binary = qd_df$Dam_binary)

# Also, need to do this for each of the 95% CIs, but the order of the
# de-scaling matters, so doing this twice more with each of the
# intervals as the first column.

# 2.5% lower interval
qd_select25 <- qd_df %>%
  dplyr::select(`lower__`, NHD_RdDensWs, log_width) %>%
  rename("lower_QcQ2" = "lower__")

qd_descaled_data25 <- as.data.frame(t(t(qd_select25) * scale3 + center3)) %>%
  mutate(Dam_binary = qd_df$Dam_binary) %>%
  dplyr::select(Dam_binary, lower_QcQ2)

# 97.5% lower interval
qd_select975 <- qd_df %>%
  dplyr::select(`upper__`, NHD_RdDensWs, log_width) %>%
  rename("upper_QcQ2" = "upper__")

qd_descaled_data975 <- as.data.frame(t(t(qd_select975) * scale3 + center3)) %>%
  mutate(Dam_binary = qd_df$Dam_binary) %>%
  dplyr::select(Dam_binary, upper_QcQ2)

qd_descaled_data <- left_join(qd_descaled_data, qd_descaled_data25)
qd_descaled_data <- left_join(qd_descaled_data, qd_descaled_data975)

(plot_qd <- ggplot(qd_descaled_data, aes(x = Dam_binary, y = 10^logQcQ2)) +
    # geom_point(size = 3) +
    # geom_errorbar(aes(ymin = lowerQcQ2, ymax = upperQcQ2), 
    #               width = 0.2) +
    geom_jitter(data = dat_Qc_trim %>%
                  drop_na(Dam_binary), aes(x = Dam_binary, y = Qc_Q2yr),
                alpha = 0.4, width = 0.1, color = "#5A7ECB") +
    labs(x = "Likelihood of Interference by Dams",
         y = expression(Q[c]:Q[2~yr])) +
    scale_x_discrete(labels = c("5-50%", "100%")) +
    scale_y_log10() +
    theme_bw())

####### Size #######

q_w <- conditional_effects(q1, effects = "log_width")

# Create new dataframe
qw_df <- q_w$log_width

qw_select <- qw_df %>%
  dplyr::select(`estimate__`, NHD_RdDensWs, `effect1__`) %>%
  rename("logQcQ2" = "estimate__",
         "log_width" = "effect1__")

# And calculate true yield values
qw_descaled_data <- as.data.frame(t(t(qw_select) * scale3 + center3))

# Also, need to do this for each of the 95% CIs, but the order of the
# de-scaling matters, so doing this twice more with each of the
# intervals as the first column.

# 2.5% lower interval
qw_select25 <- qw_df %>%
  dplyr::select(`lower__`, NHD_RdDensWs, `effect1__`) %>%
  rename("lower_QcQ2" = "lower__",
         "log_width" = "effect1__")

qw_descaled_data25 <- as.data.frame(t(t(qw_select25) * scale3 + center3)) %>%
  dplyr::select(log_width, lower_QcQ2)

# 97.5% lower interval
qw_select975 <- qw_df %>%
  dplyr::select(`upper__`, NHD_RdDensWs, `effect1__`) %>%
  rename("upper_QcQ2" = "upper__",
         "log_width" = "effect1__")

qw_descaled_data975 <- as.data.frame(t(t(qw_select975) * scale3 + center3)) %>%
  dplyr::select(log_width, upper_QcQ2)

qw_descaled_data <- left_join(qw_descaled_data, qw_descaled_data25)
qw_descaled_data <- left_join(qw_descaled_data, qw_descaled_data975)

(plot_qw <- ggplot(qw_descaled_data, aes(x = 10^log_width, y = 10^logQcQ2)) +
    # geom_line(color = "black", linewidth = 1) +
    # geom_ribbon(aes(ymin = lowerQcQ2, ymax = upperQcQ2),
    #             alpha = 0.25) +
    geom_point(data = dat_Qc_trim, aes(x = 10^log_width, y = Qc_Q2yr),
               alpha = 0.4, color = "#D46F10") +
    scale_y_log10() +
    scale_x_log10()+
    labs(x = "River Width (m)",
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

####### Combined #######

# Now, let's combine the above using patchwork.
(fig_cond_Qc <- (q_fig + plot_qd + plot_qr + plot_qw) +
    plot_layout(nrow = 2) +
    plot_annotation(tag_levels = 'A'))

# And export.
# ggsave(fig_cond_Qc,
#        filename = "figures/teton_fall22/brms_Qc_cond_033123.jpg",
#        width = 28,
#        height = 20,
#        units = "cm")

# Save out this figure.
# ggsave(q_fig,
#        filename = "figures/teton_fall22/brms_QcQ2_033123.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

# End of script.
