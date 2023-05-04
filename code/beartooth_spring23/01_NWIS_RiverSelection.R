## Resilience of Stream Productivity to Disturbance
## Originally created: October 10, 2022
## Heili Lowman

#### READ ME ####

# The following set of scripts will walk through the steps necessary to
# prep and send data to Beartooth as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_XX" folders have been ignored using git.ignore,
# so links to the raw data sets are provided below.

# If you are accessing the code via GitHub, these will need to be downloaded 
# and added to a folder of the appropriate name prior to running the code.

#### Setup ####

## Load packages
lapply(c("tidyverse", "cowplot","lubridate",
         "data.table","patchwork", "here"), require, character.only=T)

## Load datasets
## Import site metadata from Appling et al.
# https://www.sciencebase.gov/catalog/item/59bff64be4b091459a5e098b
site <- fread("data_raw/site_data.tsv")

## Secondary stream order source from hypoxia metadata set
# https://www.sciencebase.gov/catalog/item/606f60afd34ef99870188ee5
## and subsetted to Appling (PC)
hyp <- fread("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
hyp <- hyp[which(hyp$DB_Source == "PC"), c("SiteID","ORD_STRA","NHD_STREAMORDE")]
colnames(hyp)[which(colnames(hyp) == "SiteID")] <- "site_name"

## Merge hypoxia metadataset with site metadata
mdf <- left_join(site, hyp, "site_name")

## Subset for site information regarding flow impediments, order, etc.
sub_mdf <- mdf[,c("site_name","long_name",
             "site_type","struct.canal_flag","struct.dam_flag",
             "struct.npdes_flag","ORD_STRA","NHD_STREAMORDE")]

#### First Set of Model Diagnostic Filters ####

# Import and subset model diagnostics from Appling et al.
# for the sites in the above subset of sites.
# https://www.sciencebase.gov/catalog/item/59eb9bafe4b0026a55ffe382 
diagnostics <- read.table("data_raw/diagnostics.tsv",sep = "\t", header=T)
diagnostics <- diagnostics[which(diagnostics$site %in% sub_mdf$site_name),]

# Filtering for Rhat (Gelman-Rubin statistic) below 1.05, 
# because this provides an indication for how well the 
# streamMetabolizer output performed

# 433 sites at the start

# quick histograms to visualize the thresholds at which I'm filtering
hist(diagnostics$err_obs_iid_sigma_Rhat)
hist(diagnostics$err_proc_iid_sigma_Rhat)
hist(diagnostics$neg_GPP)
hist(diagnostics$pos_ER) # this looks worst

highq_sites <- diagnostics[which(
  # the following two filters remove sites at which the error terms 
  # poorly converged
                diagnostics$err_obs_iid_sigma_Rhat < 1.05 & # 101 sites drop
                diagnostics$err_proc_iid_sigma_Rhat < 1.05 & # 2 sites drop
                                     
  # the following two filters refer to "percent of GPP estimates < -0.5 
  # and the percent of ER estimates > 0.5 
  # (both are biologically unrealistic outcomes)"
  # at neg_GPP < 5/15/25/50, 58/12/7/4 more sites disappear
                diagnostics$neg_GPP < 15 &
                                     
  # at pos_ER < 5/15/25/50, 85/42/32/13 more sites disappear
                diagnostics$pos_ER < 15),]

# 276 records remaining after diagnostic filtering
# but some sites have multiple time resolutions (e.g., 15 and 30 minutes)
highq_site_names <- unique(highq_sites$site) ## 248 unique sites

# Create subset s based on high quality sites
s <- sub_mdf %>%
  filter(site_name %in% highq_site_names)

# Removing all except stream sites. Possible categories include:
# ST (stream), ST-CA (canal), ST-DCH (ditch), ST-TS (tidal stream), or SP (spring)
s <- s[which(s$site_type == "ST"),] ## 8 sites drop off

# Import time series of streamMetabolizer-generated predictions
# to be able to filter by daily Rhat values for K600
# https://www.sciencebase.gov/catalog/item/59eb9c0ae4b0026a55ffe389
NWIS <- read.table("data_raw/daily_predictions.tsv", sep='\t', header = TRUE)
NWIS$date <- ymd(NWIS$date)

# Due to date issues below, forcing the dataframe to present dates in ascending order
NWIS <- NWIS %>%
  group_by(site_name) %>%
  arrange(date) %>% # ascending is the default
  ungroup()

# Instead of filtering by site-level diagnostics above, now filtering by 
# daily Rhat values
NWIS_ed <- NWIS %>%
  filter(GPP.Rhat < 1.05) %>% # remove days on which Rhat > 1.05
  filter(K600.Rhat < 1.05) # remove days on which Rhat > 1.05

# Removed 5% of days (470,763 records remaining of the original 490,907)

## Subset columns and sites
NWIS_sub <- NWIS_ed[,c("site_name","date","GPP","GPP.lower","GPP.upper", 
                    "GPP.Rhat","ER","ER.lower","ER.upper","K600",
                    "K600.lower","K600.upper","temp.water",
                    "discharge","shortwave","velocity")]

# Rename columns
colnames(NWIS_sub) <- c("site_name","date","GPP","GPP.lower","GPP.upper", 
                        "GPP.Rhat","ER","ER.lower","ER.upper","K600",
                        "K600.lower","K600.upper","temp","Q","light",
                        "velocity")

## Subset daily diagnostic filtered sites (NWIS_sub) from site-level diagnostic
## filtered sites (s)
NWIS_sub <- NWIS_sub[which(NWIS_sub$site_name %in% s$site_name),]

# So, with site & daily diagnostic filters revised, removed 31% of data records
# 341,567 of original 490,907 daily records remaining

# Confirm remaining number of sites
length(unique(NWIS_sub$site_name)) ## 240 sites

#### Second Set of Temporal Filters ####

# The secondary set of filters focus on time-series evaluation for gaps and 
# quantity of days available.

## Identify which sites have the most continuous data
NWIS_sub$doy <- yday(NWIS_sub$date)
NWIS_sub$year <- year(NWIS_sub$date)

## count days per year
dat_per_year <- NWIS_sub %>%
  dplyr::count(site_name, year)

hist(dat_per_year$n) # most records appear fairly complete
median(dat_per_year$n) # 248.5

## identify the max day gap PER YEAR within a given site
gap_per_year <- NWIS_sub %>%
  group_by(site_name, year) %>%
  mutate(gap = doy - lag(doy, default=doy[1], order_by = doy)) %>%
  ungroup()
# The code above really struggles between the dplyr/plyr packages it seems.
# So, ALWAYS double check the output above to be sure there aren't any gaps
# that are negative, as that indicates the data isn't being grouped properly.

maxgap <- gap_per_year %>%
  group_by(site_name, year) %>%
  summarize_at(.vars = "gap", .funs = max) %>%
  ungroup()

# Also, added code below because the shortened time frames (3 months)
# need to be included prior to filtering out for long gaps

# So, the longest time gap between days may be from fall of one year
# to summer of the next, but for the sites at which we only have annual
# summer data, this is a gap we must accept. More information will be
# gained from these kinds of sites than from sites with 90days of data
# but enormous gaps within the year (i.e., 30 days in March, July, and Oct).
# That's also what our re-initialization of the model will address.

## merge with number of days per year
sub_by_gap <- merge(maxgap, dat_per_year, by=c("site_name","year")) 
# 1538 site-years

# Changing to allow for season-long data to increase the number of sites.
# Filter for at least 90 days per year (minimum possible 3 month period).
sub_by_gap1 <- sub_by_gap[which(sub_by_gap$n > 90),] # 1271 site-years
sub_by_gap_sum <- sub_by_gap1 %>% 
  group_by(site_name) %>% 
  count() %>%
  ungroup() # 231 sites

## subset for sites with a max gap of 14 days
sub_by_gap2 <- sub_by_gap1[which(sub_by_gap1$gap <= 14),] # 755 site-years
length(levels(as.factor(sub_by_gap2$site_name))) # 198 unique sites

# Renaming dataset with sites, number of days, and gap records by year
highq_years <- sub_by_gap2

## Subset dataset filtered by model diagnostics by dataset filtered by
# temporal gaps
TS <- NWIS_sub[which(NWIS_sub$site_name %in% highq_years$site_name),]
# 318,480 observations remaining

## Subset also to YEARS that meet criteria
highq_years$site_year <- paste(highq_years$site_name,highq_years$year,sep = "_")
TS$site_year <- paste(TS$site_name, TS$year,sep = "_")
TS <- TS[which(TS$site_year %in% highq_years$site_year),] # Site data
# 213,519 observations remaining

# And create an accompanying, filtered metadata file 
# from the original metadata (s)
TS_site <- s[which(s$site_name %in% TS$site_name),] # Site info

## Create a new column duplicating the median GPP
TS$GPP_temp <- TS$GPP

# Where GPP is negative, replace values with a small, randomly generated number
set.seed(148)
TS[which(TS$GPP < 0),]$GPP_temp <- sample(exp(-6):exp(-4), 1)

# Calculate site-level mean and max daily GPP
TS_gpp <- TS %>%
  group_by(site_name) %>%
  summarise_at(.vars = "GPP_temp", .funs = c(mean, max)) %>%
  ungroup()

colnames(TS_gpp) <- c("site_name","GPP_mean","GPP_max")

# And append these to the larger metadataset
TS_site <- left_join(TS_site, TS_gpp, by="site_name")

#### Trim Dataframe ####

# There are 198 unique sites in the TS_site metadataset, with these new filters.
# Joining the full names to the TS dataset to better ID them when plotting.
name_bridge <- TS_site %>%
  select(site_name, long_name)

TS <- TS %>%
  left_join(name_bridge)

# subset necessary data
# GPP data
site_subset <- TS # double-check, 213,519 observations still

# site info/metadata
TS_site_subset <- mdf[which(mdf$site_name %in% site_subset$site_name),]

# data gaps by site-year
site_subset_numdays <- rbind(highq_years[which(highq_years$site_name %in% 
                                                highq_years$site_name),])

colnames(site_subset_numdays) <- c("site_name","year","max_gap","Ndays","site_year")

#### Check Covariate Data Quality ####

site_sub_list <- split(site_subset, site_subset$site_name)

plotting_covar <- function(x) {
  
  df <- x
  p <- plot_grid(

    ggplot(df, aes(date, GPP))+
      geom_point(color="chartreuse4", size=2)+
      geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="darkolivegreen4")+
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=df$long_name[1],
           subtitle=df$site_name[1])+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12)),
    
    ggplot(df, aes(date, Q))+
      geom_line(size=1.5, color="deepskyblue4")+
      labs(y="Q (cms)")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12)),
    
    ggplot(df, aes(date, temp))+
      geom_line(size=1.5, color="#A11F22")+
      labs(y="Water Temp (C)")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12)),
    
    ggplot(df, aes(date, light))+
      geom_point(size=2, color="darkgoldenrod3")+
      labs(y="Incoming Light", x="Date")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.text = element_text(size=12),
            axis.title = element_text(size=12)),
    ncol=1, align="hv")
  
  return(p)
  
}

# plot and save out all covariate data for 198 high quality sites :)

lapply(site_sub_list, function(x) ggsave(plot = plotting_covar(x),
       filename = paste("figures/beartooth_spring23/site_covariate_plots/",
                                x$site_name[1],"covar.jpg",sep = ""), 
                                         width = 10, 
                                         height = 8))

#### Export ####

## NWIS site subset
saveRDS(site_subset, 
        "data_working/NWIS_198sites_050423.rds") # GPP data
saveRDS(TS_site_subset, 
        "data_working/NWIS_198sitesinfo_050423.rds") # Site metadata
saveRDS(site_subset_numdays,
        "data_working/NWIS_198sitesNdays_050423.rds") # Gaps in data

# NOTE: THERE IS NO DAM FILTER.

# Additional information from the site_data.xml file:
# "4. The fields struct.canal_flag, struct.dam_flag, struct.npdes_flag are flags for likely site suitability for metabolism modeling, based on site proximity to infrastructure that could affect metabolism estimates.
# 4a. The flagging focused on three feature types:  permitted National Pollution Discharge Elimination System point sources (see USEPA_NPDES in Source Citations), dams identified in the National Inventory of Dams (USACE_NID), and canals and ditches identified in the National Hydrography Database (NHDPlusV2). National data layers for these features of interest were clipped to watersheds upstream of the sites in our dataset. The geodetic distance between each stream monitoring site and the nearest upstream feature of each type was calculated using the function GenerateNearTable in the arcpy library in Python 3.6.
# 4b. The mean 80% oxygen turnover distance on each day (DO.tdist80) was pulled from "Daily metabolism estimates and predictors", and the 50th, 80th, and 95th percentiles of the daily 80% turnover distances were then computed for each site.
# 4c. Distances between the site and upstream features were compared with the percentiles of 80% O2 turnover distance. Site flags are numeric to indicate whether the nearest feature of each type was farther than the 95th, 80th, 50th, or 0th percentile; those four cases are represented by site flag values of 95, 80, 50, or 0, respectively, such that a value of 95 indicates the least probable interference from a structure of a given type."

# End of script.
