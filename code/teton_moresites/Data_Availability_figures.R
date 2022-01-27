# Figures for depicting data availability
# October 28, 2021
# Heili Lowman

# This script will be use for plotting the dataset based 
# on filtering performed in the teton_moresites/NWIS_RiverSelection.R script.


# Setup -------------------------------------------------------------------

# Load packages.
library(tidyverse)
library(lubridate)
library(data.table)
library(calecopal)
library(patchwork)
library(sf)
library(viridis)
library(maps)
library(mapproj)
library(here)

# Double check working directory
here()

# Load datasets.
# Main dataset containing site information from Appling + Koenig
sitesjoin <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")

# Dataset with data itself
site_subset <- readRDS("data_working/NWIS_207sites_subset.rds")

# Figures -----------------------------------------------------------------

# Distribution of available sites by stream order

main1 <- sitesjoin %>%
  distinct(site_name, .keep_all = TRUE) %>% # filter by unique site locations
  select(site_name, NHD_STREAMORDE) %>%
  mutate(order = fct_explicit_na(factor(NHD_STREAMORDE))) # ensure NA is a level

(storder <- ggplot(main1) + # base plot
  geom_histogram(aes(x = order, fill = order), stat = "count") + # stream order histogram
  scale_fill_manual(values = cal_palette("sbchannel", n = 10, type = "continuous")) + # custom colors
  labs(x = "Stream Order",
       y = "Site Count") +
  theme_classic() + # remove grid
  theme(legend.position = "none")) # remove legend

# Takeaway - Stream Order of 5 is now most represented.

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
(sitemap <- ggplot(states_sf) + # base plot
  geom_polygon(aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") + # map of states
  geom_point(data = sites_sf %>% 
               filter(site_name != "nwis_15298040"), # removing 1 alaska site for now
             aes(x = lon, y = lat, color = order), 
             size = 4, alpha = 0.8) + # map of sites
  scale_color_manual(name = "Stream Order", 
                     values = cal_palette("sbchannel", n = 10, type = "continuous")) + # custom colors
  theme_classic() + # remove grid
  labs(x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "bottom") + # reposition legend
  coord_map(projection = "albers", lat0 = 39, lat1 = 45))   

# Takeaways - Still skewed towards northeastern/midwestern US.

# Distribution of available sites by years of data

# Create intermediate dataset to summarize by year
# Note, data has already been filtered, so for the years included, data is AOK
main3 <- site_subset %>%
  mutate(yearf = factor(year)) %>%
  group_by(site_name, yearf) %>%
  summarize(meanGPP = mean(GPP)) %>%
  ungroup()

# Export for shiny app use:
saveRDS(main3, "data_working/teton_207rivers_sitesyrsgpp.rds")

# Count number of years (observations) per site
main4 <- main3 %>%
  count(site_name) %>%
  mutate(n_f = factor(n))

# Export for rmarkdown use:
saveRDS(main4, "data_working/teton_207rivers_sitesyears.rds")

(styears <- ggplot(main4) + # base plot
  geom_histogram(aes(x = n_f, fill = n_f), stat = "count") + # years data histogram
  scale_fill_manual(values = cal_palette("desert", n = 9, type = "continuous")) + # custom colors
  labs(x = "Years of Available Data",
       y = "Site Count") +
  theme_classic() + # remove grid
  theme(legend.position = "none")) # remove legend

# Combine all figures into one

full_fig <- (storder | styears) /
  sitemap  +
  plot_annotation(tag_levels = 'A')

full_fig

# ggsave(("figures/teton_moresites/map_fig_207sites.png"),
#        width = 16,
#        height = 20,
#        units = "cm"
# )

# Time Sequences ------------------------------------------------------------

# Borrowing some code from Joanna to calculate sequences of dates

# create placeholder columns for day to day differences and sequences
df <- site_subset %>%
  mutate(diff_time = 0, seq = 1)

# Need to add some more groupings since I'm iterating over more than 1 site...

# First, I'm going to get this working on a single site
d1 <- df %>%
  filter(site_name == "nwis_01124000")

d <- d1

# set the first day to 0
#d$diff_time[1] <- 0

# calculate the difference from one day to the next
for(i in 2:nrow(d)){
  d$diff_time[i] = difftime(time1 = d$date[i], time2 = d$date[(i-1)],
                            units = "days")
}

# convert to character
#d$diff_time <- as.character(as.numeric(d$diff_time))

# set the first in the sequence to 1
#d$seq[1] <- 1

# delineate sequenced time frames based on day to day differences
# anything less than 14 day gaps is permissible
for(i in 2:nrow(d)){
  if(d$diff_time[i] < 14){
    d$seq[i] = d$seq[(i-1)]
  } else {
    d$seq[i] = d$seq[(i-1)] + 1
  }
}

# split into a list based on events
l1 <- split(d, as.factor(d$seq))

# create larger dataframe as list
my_list <- split(df, f = df$site_name)

# create function for application to all sites
seqFUN <- function(d){
  
  # calculate the difference from one day to the next
  for(i in 2:nrow(d)){
    d$diff_time[i] = difftime(time1 = d$date[i], time2 = d$date[(i-1)],
                              units = "days")
  }
  
  # delineate sequenced time frames based on day to day differences
  # anything less than 14 day gaps is permissible
  for(i in 2:nrow(d)){
    if(d$diff_time[i] < 14){
      d$seq[i] = d$seq[(i-1)]
    } else {
      d$seq[i] = d$seq[(i-1)] + 1
    }
  }
  
  # split into a list based on events
  l <- split(d, as.factor(d$seq))
  
  return(l)
  
}

# And now map this to the entire site list.
events_dat <- lapply(my_list, seqFUN)

# trying out a first figure to try and visualize individual events
# first, make the output list back into a dataframe
# need to unnest the doubly nested list, so this is the first round
events_dat1 <- lapply(events_dat, rbindlist)

# and then unnest the final list level
events_df <- rbindlist(events_dat1)

# calculate some summary stats for plotting
events_summary <- events_df %>%
  count(site_name, seq)

mean(events_summary$n) #458

(evyears <- ggplot(events_summary) + # base plot
    geom_point(aes(x = site_name, y = n, color = site_name), size = 2) +
    scale_color_manual(values = cal_palette("collinsia", n = 207, type = "continuous")) +
    labs(x = "Sites",
         y = "Event Lengths (days)") +
    geom_hline(yintercept = 458, linetype = "dotted") +
    geom_hline(yintercept = 365, linetype = "solid") +
    theme_classic() + # remove grid
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")) # remove legend

# ggsave(("figures/teton_moresites/events_fig_207sites.png"),
#        width = 20,
#        height = 20,
#        units = "cm"
# )

# some more summary stats for plotting
events_summary2 <- events_summary %>%
  group_by(site_name) %>%
  summarize(max_seq = max(seq, na.rm=TRUE), mean_event = mean(n, na.rm=TRUE)) %>%
  ungroup()

(evcounts <- ggplot(events_summary2) + # base plot
    geom_point(aes(x = max_seq, y = mean_event, color = site_name), size = 2) +
    scale_color_manual(values = cal_palette("collinsia", n = 207, type = "continuous")) +
    labs(x = "Max. # of Events",
         y = "Mean Event Length (days)") +
    geom_hline(yintercept = 458, linetype = "dotted") +
    geom_hline(yintercept = 365, linetype = "solid") +
    theme_classic() + # remove grid
    theme(#axis.text.x = element_blank(),
          #axis.ticks.x = element_blank(),
          legend.position = "none")) # remove legend

# ggsave(("figures/teton_moresites/events2_fig_207sites.png"),
#        width = 20,
#        height = 20,
#        units = "cm"
# )

# Creating an additional figure for use in my presentation to Modelscape.
# Import time series of streamMetabolizer-generated predictions
NWIS <- read.table("data_raw/daily_predictions.tsv", sep='\t', header = TRUE)

wv <- NWIS %>%
  filter(site_name == "nwis_01608500")# filter only for stream in WV

wv$date <- ymd(wv$date) # and structure dates properly

wv2012 <- wv %>%
  mutate(year = year(date)) %>%
  filter(year ==2012)

wv12_do <- ggplot(wv2012, aes(date, DO.obs))+
  geom_line(color="gray60", size=1)+
  labs(y=expression('DO (mg '*~O[2]~ L^-1*')'), x = "Date")+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12))
  
wv12_gpp <- ggplot(wv2012, aes(date, GPP))+
  geom_point(color="chartreuse4", size=2)+
  geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="darkolivegreen4")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x = "Date")+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12))

wv12_q <- ggplot(wv2012, aes(date, discharge))+
  geom_line(size=1.5, color="deepskyblue4")+
  labs(y=expression('Q (cm '*s^-1*')'), x = "Date")+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12))

wv12fig <- wv12_do / wv12_gpp / wv12_q

# ggsave(("figures/presentations/nwis_01608500_2012_do_gpp_q.png"),
#        width = 20,
#        height = 20,
#        units = "cm"
# )

# Also, creating a figure showing time gaps to visualize reinitialization.

pa <- NWIS %>%
  filter(site_name == "nwis_03007800")# filter only for stream in PA

pa$date <- ymd(pa$date) # and structure dates properly

pa <- pa %>%
  mutate(year = year(date)) %>%
  filter(year < 2015)

pa_gpp <- ggplot(pa, aes(date, GPP))+
  geom_point(color="chartreuse4", size=2)+
  geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="darkolivegreen4")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x = "Date")+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12))

# ggsave(("figures/presentations/nwis_03007800_2011_2014_gpp.png"),
#        width = 25,
#        height = 8,
#        units = "cm"
# )

# End of script.
