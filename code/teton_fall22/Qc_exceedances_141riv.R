## Estimating discharge disturbance threshold exceedances
## December 1, 2022
## Heili Lowman

# The following script will calculate the number of times each site
# exceeds what we have estimated as the critical disturbance threshold (Qc).

# Then, it will examine if the number of exceedances has any relationship
# estimated rmax values.

# Load packages.
lapply(c("lubridate","tidyverse", "here", "viridis",
         "reshape2","ggExtra","patchwork"), require, character.only=T)

# Load datasets.

# First, the raw data loaded into the model.
dat_in <- readRDS("data_working/list_182sites_Qmaxnorm_allSL.rds")

# And then a dataset containing Qc values that have been converted from
# c values estimated by our model.
dat_Qc <- readRDS("data_working/QcQ2_159sites_120822.rds")

# Finally the dataset containing rmax & amax values for plotting purposes.
dat_rmax <- readRDS("data_working/rmax_filtered_159sites_113022.rds")
dat_amax <- readRDS("data_working/maxalgalyield_159sites_021323.rds")

# Also also, hypoxia dataset for additional info re: slope
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")

# Now, to create a practice workflow to later turn into a function to
# quantify exceedances of Qc. Going to first pick a site where it actually
# happens.

test_df <- dat_in$nwis_01400500

# First visualize the data to see how many times it happens
ggplot(test_df, aes(x = Date, y = Q)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 67.4526324) +
  labs(x = "Date",
       y = "Discharge") +
  theme_bw()
# Ok, so roughly 20 times.

# For ease, making the original list into a df.
dat_in_df <- map_df(dat_in, ~as.data.frame(.x), .id="site_name")

# Trim down to columns of interest for Qc dataset
dat_Qc_trim <- dat_Qc %>%
  dplyr::select(site_name, Qc)

# Need to add Qc values for each site to the input data
dat_in_Qc <- inner_join(dat_in_df, dat_Qc_trim)

# And make back into a list to apply the function below.
## split list by ID
dat_in_Qc_l <- split(dat_in_Qc, dat_in_Qc$site_name)

# Modifying my event delineation function from step3 of this workflow.

# Event delineation function:
# loop over the separate time sequences for a given site
exceedFUN <- function(d){
  
  # calculate the difference from one day to the next - YES/NO
    
    d <- d %>%
      mutate(exceedance = case_when(Q > Qc ~ "YES", TRUE ~ "NO"))
  
  # delineate sequenced time frames based on day to day differences
    
    d <- d %>% # mark all yesses preceded by nos
      mutate(sequence = case_when(exceedance == "YES" & 
                               lag(exceedance) == "NO" ~ 1,
                             TRUE ~ 0)) %>%
      # and count total days of exceedance
      mutate(total = case_when(exceedance == "YES" ~ 1,
                               TRUE ~ 0))
  
  return(d)
  
}

# apply event delineation function
dat_in_exc <- lapply(dat_in_Qc_l, function(x) exceedFUN(x))

# spot checking a site to see how it looks
site1 <- dat_in_exc$nwis_02203900

# Make the list back into a df
dat_in_exc_df <- map_df(dat_in_exc, ~as.data.frame(.x), .id="site_name")

# And summarize total exceedance events and exceedance days as well as
# total days in a record.
dat_exceed_sum <- dat_in_exc_df %>%
  group_by(site_name) %>%
  summarize(total_exc_events = sum(sequence),
            total_exc_days = sum(total),
            total_days = n()) %>%
  ungroup() %>%
  # calculate # of exceedance days relative to length of dataset
  mutate(total_exc_rel = total_exc_days/total_days)

# Join with rmax values.
dat_exc_rmax <- inner_join(dat_exceed_sum, dat_rmax)

# And plot to examine for relationships.
(exceed_fig1 <- ggplot(dat_exc_rmax, aes(x = total_exc_events, y = r_med)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = expression(Total~Q[c]~Exceedance~Events),
       y = expression(Maximum~Growth~Rate~(r[max]))) +
  theme_bw() +
  theme(legend.position = "none"))

(exceed_fig1a <- ggplot(dat_exc_rmax, aes(x = total_exc_events, y = r_med)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(x = expression(Total~Q[c]~Exceedance~Events),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    scale_x_log10() + 
    theme_bw() +
    theme(legend.position = "none"))

(exceed_fig2 <- ggplot(dat_exc_rmax, aes(x = total_exc_days, y = r_med)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(x = expression(Total~Q[c]~Exceedance~Days),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw() +
    theme(legend.position = "none"))

(exceed_fig2a <- ggplot(dat_exc_rmax, aes(x = total_exc_days, y = r_med)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(x = expression(Total~Q[c]~Exceedance~Days),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    scale_x_log10() + 
    theme_bw() +
    theme(legend.position = "none"))

(fig_exc <- exceed_fig1 + exceed_fig1a +
    exceed_fig2 + exceed_fig2a +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 2))

# And export for use in the Rmarkdown file.
# ggsave(fig_exc,
#        filename = "figures/teton_fall22/rmax_exceedance_120822.jpg",
#        width = 20,
#        height = 20,
#        units = "cm") # n = 159

# Additional figures after speaking with Joanna.
(exceed_fig3 <- ggplot(dat_exc_rmax, aes(x = total_exc_rel*100, 
                                         y = r_med,
                                         color = width_med)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(x = expression(Percent~Q[c]~Exceedance~Days~'in'~Record),
         y = expression(Maximum~Growth~Rate~(r[max])),
         color = "Stream Width") +
    scale_color_viridis() +
    scale_x_log10() + 
    theme_bw())

# Trim site info dataset for slope info only
site_info_trim <- site_info %>%
  select(SiteID, slope_calc)

# Join with exceedance and rmax values.
dat_exc_rmax_slope <- left_join(dat_exc_rmax, site_info_trim, by = c("site_name" = "SiteID"))

# And additional figures using slope.
(exceed_fig4 <- ggplot(dat_exc_rmax_slope, aes(x = total_exc_rel*100, 
                                         y = slope_calc,
                                         color = width_med)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(x = expression(Percent~Q[c]~Exceedance~Days~'in'~Record),
         y = "Slope",
         color = "Stream Width") +
    scale_color_viridis() +
    scale_x_log10() + 
    theme_bw())

# Pulling out Qc:Q2yr values and adding to dataset.
dat_Qc_trim2 <- dat_Qc %>%
  select(site_name, Qc_Q2yr)

dat_exc_rmax_slope_QcQ2 <- left_join(dat_exc_rmax_slope, dat_Qc_trim2)

(exceed_fig4.2 <- ggplot(dat_exc_rmax_slope_QcQ2, aes(x = Qc_Q2yr, 
                                               y = slope_calc,
                                               color = width_med)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(x = expression(Q[c]~to~Q["2yr"]),
         y = "Slope",
         color = "Stream Width") +
    scale_color_viridis() +
    scale_x_log10() + 
    theme_bw())

(fig_exc2 <- exceed_fig3 + exceed_fig4.2 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 1))

# And export for use in the Rmarkdown file.
# ggsave(fig_exc2,
#        filename = "figures/teton_fall22/rmax_slope_exceedance_021423.jpg",
#        width = 30,
#        height = 10,
#        units = "cm") # n = 159

# And export just the exceedance data for use in the linear models.
saveRDS(dat_exceed_sum, "data_working/Qc_exceedances_159sites_021423.rds")

# One more figure to examine accrual vs. exceedance.

# Join with amax values.
dat_exc_amax <- inner_join(dat_exceed_sum, dat_amax)

(exceed_fig0a <- ggplot(dat_exc_amax, aes(x = total_exc_events, 
                                          y = yield_med2)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(x = expression(Total~Q[c]~Exceedance~Events),
         y = expression(Maximum~Biomass~Accrual~(a[max]))) +
    scale_x_log10() + 
    #scale_y_log10() + 
    theme_bw() +
    theme(legend.position = "none"))

# Standardize by years on record for each site.
dat_years <- dat_in_df %>%
  group_by(site_name) %>%
  summarize(years = n_distinct(year)) %>%
  ungroup()

# Join with accrual dataset.
dat_exc_amax_y <- inner_join(dat_exc_amax, dat_years)

# and make a column with exceedances normalized by # of years on record.
dat_exc_amax_y <- dat_exc_amax_y %>%
  mutate(exc_y = total_exc_events/years)

(exceed_fig0b <- ggplot(dat_exc_amax_y, aes(x = exc_y, 
                                          y = yield_med2)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(x = expression(Total~Q[c]~Exceedance~Events~per~Year),
         y = expression(Maximum~Biomass~Accrual~(a[max]))) +
    scale_x_log10() + 
    #scale_y_log10() + 
    theme_bw() +
    theme(legend.position = "none"))

# End of script.
