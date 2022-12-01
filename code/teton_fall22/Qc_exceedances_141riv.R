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
# c values estimated by our model and filtered for high-performing sites.
dat_Qc <- readRDS("data_working/QcQ2_filtered_141sites_113022.rds")

# Finally the dataset containing rmax values for plotting purposes.
dat_rmax <- readRDS("data_working/rmax_filtered_159sites_113022.rds")

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
  select(site_name, Qc)

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
      mutate(exceedance = case_when(Q > Qc ~ "YES", TRUE~ "NO"))
  
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

# And summarize total exceedance events and exceedance days
dat_exceed_sum <- dat_in_exc_df %>%
  group_by(site_name) %>%
  summarize(total_exc_events = sum(sequence),
            total_exc_days = sum(total)) %>%
  ungroup()

# Join with rmax values.
dat_exc_rmax <- inner_join(dat_exceed_sum, dat_rmax)

# And plot to examine for relationship.
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
#        filename = "figures/teton_fall22/rmax_exceedance_120122.jpg",
#        width = 20,
#        height = 20,
#        units = "cm") # n = 141

# And export just the exceedance data for use in the linear models.
saveRDS(dat_exceed_sum, "data_working/Qc_exceedances_141sites_120122.rds")

# End of script.
