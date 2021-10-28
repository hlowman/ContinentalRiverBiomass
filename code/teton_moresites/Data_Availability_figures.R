# Figures for depicting data availability
# October 28, 2021
# Heili Lowman

# This script will be use for plotting the dataset based 
# on filtering performed in the teton_moresites/NWIS_RiverSelection.R script.


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
sitesjoin <- readRDS("data_working/NWIS_73sitesinfo_subset.rds")

# Dataset with data itself
site_subset <- readRDS("data_working/NWIS_73sites_subset.rds")

# Figures -----------------------------------------------------------------

# Distribution of available sites by stream order

main1 <- sitesjoin %>%
  distinct(site_name, .keep_all = TRUE) %>% # filter by unique site locations
  select(site_name, NHD_STREAMORDE) %>%
  mutate(order = fct_explicit_na(factor(NHD_STREAMORDE))) # ensure NA is a level

storder <- ggplot(main1) + # base plot
  geom_histogram(aes(x = order, fill = order), stat = "count") + # stream order histogram
  scale_fill_manual(values = cal_palette("sbchannel", n = 9, type = "continuous")) + # custom colors
  labs(x = "Stream Order",
       y = "Site Count") +
  theme_classic() + # remove grid
  theme(legend.position = "none") # remove legend

storder

# Takeaway - Stream Order of 2 is still most represented.

# Distribution of available sites by geography

main2 <- sitesjoin %>%
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

# Takeaways - Still skewed towards northeastern US, but this was true for the past run.

# Distribution of available sites by years of data

# Create intermediate dataset to summarize by year
# Note, data has already been filtered, so for the years included, data is AOK
main3 <- site_subset %>%
  mutate(yearf = factor(year)) %>%
  group_by(site_name, yearf) %>%
  summarize(meanGPP = mean(GPP)) %>%
  ungroup()

# Export for shiny app use:
saveRDS(main3, "data_working/teton_73rivers_sitesyrsgpp.rds")

# Count number of years (observations) per site
main4 <- main3 %>%
  count(site_name) %>%
  mutate(n_f = factor(n))

# Export for rmarkdown use:
saveRDS(main4, "data_working/teton_73rivers_sitesyears.rds")

styears <- ggplot(main4) + # base plot
  geom_histogram(aes(x = n_f, fill = n_f), stat = "count") + # years data histogram
  scale_fill_manual(values = cal_palette("desert", n = 9, type = "continuous")) + # custom colors
  labs(x = "Years of Available Data",
       y = "Site Count") +
  theme_classic() + # remove grid
  theme(legend.position = "none") # remove legend

styears

# Combine all figures into one

full_fig <- (storder | styears) /
  sitemap  +
  plot_annotation(tag_levels = 'A')

full_fig

# ggsave(("figures/teton_moresites/map_fig_73sites.png"),
#        width = 16,
#        height = 20,
#        units = "cm"
# )

# End of script.
