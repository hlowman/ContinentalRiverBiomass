# Data Availability
# April 28, 2021
# Heili Lowman

# This script will be use for preliminary plotting
# of the Appling et al. 2018 dataset based on filtering
# performed in the data_qaqc_prelim.R script.

# Setup -------------------------------------------------------------------

# Load packages.
library(tidyverse)
library(lubridate)
library(calecopal)
library(patchwork)
library(maps)
library(here)

# Double check working directory
here()

# Load datasets.
# Main dataset containing site information from Appling + Koenig
sitesjoin <- readRDS(here("data_working", "sitesjoin.rds"))

# Main dataset containing filtered Appling dataset
filterdat <- readRDS(here("data_working", "filterdat.rds"))


# Tidy --------------------------------------------------------------------

# Join two main datasets together
main <- left_join(filterdat, sitesjoin, by = "site_name")


# Figures -----------------------------------------------------------------

# Distribution of available sites by stream order

main1 <- main %>%
  distinct(site_name, .keep_all = TRUE) %>% # filter by unique site locations
  select(site_name, NHD_STREAMORDE) %>%
  mutate(order = fct_explicit_na(factor(NHD_STREAMORDE))) # ensure NA is a level

storder <- ggplot(main1) + # base plot
  geom_histogram(aes(x = order, fill = order), stat = "count") + # stream order histogram
  scale_fill_manual(values = cal_palette("kelp1", n = 9, type = "continuous")) + # custom colors
  labs(x = "Stream Order",
       y = "Site Count") +
  theme_classic() + # remove grid
  theme(legend.position = "none") # remove legend

storder

# Takeaway - Stream Order of 5 is most represented, otherwise somewhat even from 1-6.
# 10 or less sites of stream order 7, 9, and not classified.

# Distribution of available sites by geography

main2 <- main %>%
  distinct(site_name, .keep_all = TRUE) %>% # filter by unique site locations
  mutate(site_typef = factor(site_type, levels = c("SP", "ST-TS", "ST", "ST-CA", "ST-DCH")),
         order = fct_explicit_na(factor(NHD_STREAMORDE))) # adding levels by site type and order

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

# site map colored by site type
sitemap2 <- ggplot(states_sf) + # base plot
  geom_polygon(aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") + # map of states
  geom_point(data = sites_sf, aes(x = lon, y = lat, color = site_typef), 
             size = 2, alpha = 0.75) + # map of sites
  scale_color_manual(name = "Site Type", 
                     values = c("royalblue4", "steelblue", 
                                "lightskyblue2", "lightpink", "red2")) + # custom colors
  theme_classic() # remove grid

sitemap2

# site map colored by stream order
sitemap3 <- ggplot(states_sf) + # base plot
  geom_polygon(aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") + # map of states
  geom_point(data = sites_sf, aes(x = lon, y = lat, color = order), 
             size = 2, alpha = 0.75) + # map of sites
  scale_color_manual(name = "Stream Order", 
                     values = cal_palette("kelp1", n = 9, type = "continuous")) + # custom colors
  theme_classic() # remove grid

sitemap3

# Takeaways - Still skewed towards eastern US, but this was true for the overall dataset. Data now only available for contiguous US. Stream orders appear to be distributed fairly evenly.

# Distribution of available sites by years of data

# Create intermediate dataset to summarize by year
main3 <- main %>%
  mutate(yearf = factor(year)) %>%
  group_by(site_name, yearf) %>%
  summarize(meanGPP = mean(GPP)) %>%
  ungroup()

# Count number of years (observations) per site
main4 <- main3 %>%
  count(site_name) %>%
  mutate(n_f = factor(n))

styears <- ggplot(main4) + # base plot
  geom_histogram(aes(x = n_f, fill = n_f), stat = "count") + # years data histogram
  scale_fill_manual(values = cal_palette("sbchannel", n = 9, type = "continuous")) + # custom colors
  labs(x = "Years of Available Data",
       y = "Site Count") +
  theme_classic() + # remove grid
  theme(legend.position = "none") # remove legend

styears

# Combine all figures into one

full_fig <- (storder | styears) /
  sitemap2 /
  sitemap3

full_fig + plot_annotation(title = "Data Availability After Initial Filtering")

# ggsave(("figures/full_fig.png"),
#        width = 15,
#        height = 20,
#        units = "cm"
# )

# End of script.
