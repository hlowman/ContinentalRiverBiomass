# Data Exploration
# April 23, 2021
# Heili Lowman

# This script will be use for *very* preliminary exploration
# of the Appling et al. 2018 dataset (citation below).

# Appling, A.P., Read, J.S., Winslow, L.A., Arroita, M.,
# Bernhardt, E.S., Griffiths, N.A., Hall, R.O., Jr., Harvey, 
# J.W., Heffernan, J.B., Stanley, E.H., Stets, E.G., and
# Yackulic, C.B., 2018, Metabolism estimates for 356 U.S.
# rivers (2007-2017): U.S. Geological Survey data release,
# https://doi.org/10.5066/F70864KX.


# Setup -------------------------------------------------------------------

# Load packages.
library(tidyverse)
library(calecopal)
library(sf)
library(maps)
library(mapproj)
library(lubridate)
library(viridis)
library(here)

# Double check working directory
here()

# Load datasets.
sitedat <- read_tsv(here("data_raw", 
                         "site_data.tsv"))
dailydat <- read_tsv(here("data_raw", 
                     "daily_predictions.tsv"))


# Site Data ---------------------------------------------------------------

# Initial observations:
# Site data - includes info re: canals, dams, NPDES discharge sites on a 0-100 scale
# In addition to COMID, lat/long, altitude data.
# There also appears to be multiple site types.

# Create new factored site type column to investigate possible site types.
sitedat1 <- sitedat %>%
  mutate(site_typef = factor(site_type, levels = c("SP", "ST-TS", "ST", "ST-CA", "ST-DCH"))) # adding levels in once groupings were determined

# Examine what types there are.
levels(sitedat1$site_typef)
# per https://help.waterdata.usgs.gov/site_tp_cd
# SP - spring (where water table intersects the land surface)
# ST - stream (entirely natural or partially engineered)
# ST-CA - canal (artificial, connecting two bodies of water)
# ST-DCH - ditch (artificial, smaller than a canal)
# ST-TS - tidal stream (flow is influenced by tide, but water chemistry is not)

# Will we be focused only on ST? Or perhaps separate out results by natural/artificial streams? Per meeting 4/23/2021, not worried about filtering out any sites based on site type for now.

# There's a WIDE altitudinal gradient
min(sitedat$alt, na.rm = TRUE) # -14.36 feet
max(sitedat$alt, na.rm = TRUE) # 7,380 feet
# May need to convert this to meters.

# Create map of sites.

# make data sf object
site_sf <- st_as_sf(sitedat1,
                    coords = c("lon", "lat"), # always put lon (x) first
                    remove = F, # leave the lat/lon columns in too
                    crs = 4269) # projection: NAD83

# make base US map using code from https://www.nceas.ucsb.edu/sites/default/files/2020-04/OverviewCoordinateReferenceSystems.pdf
states <- map_data("state")
state_sf <- st_as_sf(states,
                     coords = c("long", "lat"),
                     remove = F,
                     crs = 4269)

sitemap <- ggplot(state_sf) + # base plot
  geom_polygon(aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") + # map of states
  geom_point(data = site_sf, aes(x = lon, y = lat, color = site_typef), 
             size = 2, alpha = 0.75) + # map of sites
  scale_color_manual(name = "Site Type", 
                     values = c("royalblue4", "steelblue", 
                                "lightskyblue2", "lightpink", "red2")) + # custom colors
  theme_classic() # remove grid

sitemap

# Additional insights from mapping:
# All states represented except AZ, HI, NH, NV, VT
# Two sites in AK and one in PR
# Heavily skewed towards east coast/midwest sites
# Nearly all data are streams, not 100% engineered channels or springs

# Examine sites by stream order (data assembled by Lauren Koenig)
# data
# https://drive.google.com/file/d/1zMhad9X-QmdG1IC0qjMD9Tk3E1PPTAJZ/view
# metadata
# https://drive.google.com/file/d/1-FEEfTOByMPWBdeGuLiZFQadmO07-Aym/view

# Prediction Data ---------------------------------------------------------

# To first bite off a smaller chunk and examine the daily predictions,
# filter by Brandywine Creek NWIS 01481500
brandydat <- dailydat %>%
  filter(site_name == "nwis_01481500")

# And generate a quick plot to see how things look.
plot(GPP ~ date, data = brandydat, main = "Brandywine Creek GPP")

# filter + plot Sligo Creek data
sligodat <- dailydat %>%
  filter(site_name == "nwis_01650800")
plot(discharge ~ date, data = sligodat, main = "Sligo Creek Discharge")

# how many sites are there?
sitedat2 <- sitedat1 %>%
  mutate(site_namef = factor(site_name))
levels(sitedat2$site_namef) # 602 with data

dailydat1 <- dailydat %>%
  mutate(site_namef = factor(site_name))
levels(dailydat1$site_namef) # 356 with model results

# Quick and dirty temporal check
dailydat2 <- dailydat1 %>%
  mutate(datez = ymd(date)) %>%
  mutate(year = year(datez)) # format dates using lubridate

annualdat <- dailydat2 %>%
  group_by(site_namef, year) %>%
  summarize(meanGPP = mean(GPP, na.rm = TRUE)) %>%
  ungroup() # calculate mean GPP by site & year

coverage_plot <- ggplot(annualdat, aes(x = year, y = site_namef)) + # base plot
  geom_line(aes(color = site_namef)) + # add line for every year for which we have data at a given site
  labs(x = "Year",
       y = "Site") + # label axes
  scale_color_viridis(discrete=TRUE) + # custom color scale
  theme_bw() + # remove grid
  theme(axis.text.y = element_blank(), legend.position = "none") # remove site names and legend

coverage_plot

# 2011 really seems to be a turning point in terms of data availability.
# Roughly half of the sites should have at least 5 years worth of data.
# And roughly 100 sites appear to span 10 years (although that's nothing to say re: data quality).

# End of script.
