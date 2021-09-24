# Data Availability
# September 9, 2021
# Heili Lowman

# This script will be use for plotting the dataset based 
# on filtering performed in the DataSource_34rivers_StreamLight.R script.


# Setup -------------------------------------------------------------------

# Load packages.
library(tidyverse)
library(lubridate)
library(calecopal)
library(patchwork)
library(sf)
library(viridis)
library(maps)
library(here)

# Double check working directory
here()

# Load datasets.
# Main dataset containing site information from Appling + Koenig
sitesjoin <- readRDS("data_working/NWIS_34sitesinfo_subset.rds")

# Main dataset containing filtered Appling dataset
filterdat <- readRDS("data_working/df_34sites.rds")


# Tidy --------------------------------------------------------------------

# Make larger list dataset a dataframe so that it can be joined with site info
filterdat1 <- map_df(filterdat, ~as.data.frame(.x), .id="site_name")

# Join two main datasets together
main <- left_join(filterdat1, sitesjoin, by = "site_name")


# Figures -----------------------------------------------------------------

# Distribution of available sites by stream order

main1 <- main %>%
  distinct(site_name, .keep_all = TRUE) %>% # filter by unique site locations
  select(site_name, NHD_STREAMORDE) %>%
  mutate(order = fct_explicit_na(factor(NHD_STREAMORDE))) # ensure NA is a level

storder <- ggplot(main1) + # base plot
  geom_histogram(aes(x = order, fill = order), stat = "count") + # stream order histogram
  scale_fill_manual(values = cal_palette("sbchannel", n = 7, type = "continuous")) + # custom colors
  scale_y_continuous(breaks=seq(0, 10, 1)) + # fix y axis labels
  labs(x = "Stream Order",
       y = "Site Count") +
  theme_classic() + # remove grid
  theme(legend.position = "none") # remove legend

storder

# Takeaway - Stream Order of 2 is most represented.

# Distribution of available sites by geography

main2 <- main %>%
  distinct(site_name, .keep_all = TRUE) %>% # filter by unique site locations
  mutate(order = fct_explicit_na(factor(NHD_STREAMORDE))) # adding level by order

# make data sf object
sites_sf <- st_as_sf(main2,
                    coords = c("lon", "lat"), # always put lon (x) first
                    remove = F, # leave the lat/lon columns in too
                    crs = 4269) # projection: NAD83

# make base US map
states <- map_data("state")
states_sf <- st_as_sf(states,
                     coords = c("long", "lat"),
                     remove = F,
                     crs = 4269)

# site map colored by stream order
sitemap <- ggplot(states_sf) + # base plot
  geom_polygon(aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") + # map of states
  geom_point(data = sites_sf, aes(x = lon, y = lat, color = order), 
             size = 4, alpha = 0.8) + # map of sites
  scale_color_manual(name = "Stream Order", 
                     values = cal_palette("sbchannel", n = 9, type = "continuous")) + # custom colors
  theme_classic() + # remove grid
  labs(x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "bottom") # reposition legend

sitemap

# Takeaways - Still skewed towards eastern US, but this was true for the overall dataset. Data now only available for contiguous US. Stream orders appear to be distributed fairly evenly.

# Distribution of available sites by years of data

# Create intermediate dataset to summarize by year
# Note, data has already been filtered, so for the years included, data is AOK
main3 <- main %>%
  mutate(yearf = factor(year)) %>%
  group_by(site_name, yearf) %>%
  summarize(meanGPP = mean(GPP)) %>%
  ungroup()

# Export for shiny app use:
saveRDS(main3, "data_working/teton_34rivers_sitesyrsgpp.rds")

# Count number of years (observations) per site
main4 <- main3 %>%
  count(site_name) %>%
  mutate(n_f = factor(n))

# Export for rmarkdown use:
saveRDS(main4, "data_working/teton_34rivers_sitesyears.rds")

styears <- ggplot(main4) + # base plot
  geom_histogram(aes(x = n_f, fill = n_f), stat = "count") + # years data histogram
  scale_fill_manual(values = cal_palette("desert", n = 8, type = "continuous")) + # custom colors
  labs(x = "Years of Available Data",
       y = "Site Count") +
  scale_y_continuous(breaks=seq(0, 10, 1)) + # fix y axis labels
  theme_classic() + # remove grid
  theme(legend.position = "none") # remove legend

styears

# Combine all figures into one

full_fig <- (storder | styears) /
  sitemap  +
  plot_annotation(tag_levels = 'A')

full_fig

# ggsave(("figures/teton_34sites/map_fig_34sites.png"),
#        width = 16,
#        height = 20,
#        units = "cm"
# )


# Gap Investigation -------------------------------------------------------

# Additional data investigation re: time spans and gaps on June 2, 2021.

# Re-creating the coverage plot from the data_exploration script for starters.
coverage_plot <- main3 %>%
  mutate(site_namef = factor(site_name),
         year = as.numeric(as.character(yearf))) %>%
  ggplot(aes(x = year, y = site_namef)) + # base plot
  geom_line(aes(color = site_namef)) + # add line for every year for which we have data at site
  labs(x = "Year",
       y = "Site") + # label axes
  scale_color_viridis(discrete=TRUE) + # custom color scale
  theme_bw() + # remove grid
  theme(axis.text.y = element_blank(), legend.position = "none") # remove site names and legend

coverage_plot # It seems the sites with gaps are not displaying. Let's look into that more.

# Create a table of the main_4/styears histogram.
years_table <- main4 %>%
  group_by(n_f) %>%
  summarize(number_of_sites = n()) %>%
  ungroup()

# Determine which of these records are continuous.
gaps <- main3 %>%
  mutate(year = as.numeric(as.character(yearf))) %>% # recreate numeric year column
  group_by(site_name) %>%
  mutate(gap = year - lag(year, default = year[1])) %>% # calculate gaps
  ungroup() 

sites_w_gaps <- gaps %>%
  filter(gap > 1) %>% # filter only for values that denote a gap
  select(site_name) %>% # pull out only site names
  unique()

# Filter the main4 dataset (since there's one record per site) 
# to remove the sites_w_gaps determined above.
main4_cont <- anti_join(main4, sites_w_gaps)

# Create a new table.
cont_years_table <- main4_cont %>%
  group_by(n_f) %>%
  summarize(number_of_cont_sites = n()) %>%
  ungroup()

# And add to original table.
years_table <- years_table %>%
  mutate(number_of_cont_sites = cont_years_table$number_of_cont_sites)

# Takeaways:
# (1) There are 76 sites, of a total possible 202, with time gaps of some kind.
# (2) 1 year gaps are most common (n = 58 a.k.a. 75%), although there is a gap of 6 years max.
# ? Perhaps this means we should impose a filter of gaps no longer than 1 year, similar
# to how we impose the filter of gaps no longer than 2 weeks at a smaller scale.
# (3) The following table shows the number of unique sites available:
# years of data | total | continuous
# 1             | 62    | 62
# 2             | 35    | 24
# 3             | 25    | 11
# 4             | 24    | 13
# 5             | 16    | 4
# 6             | 13    | 1
# 7             | 9     | 1
# 8             | 11    | 3
# 9             | 7     | 7

# End of script.
