## Resilience of Stream Productivity to Disturbance
## July 25, 2022
## Heili Lowman

# The following set of scripts will walk through the steps necessary to
# prep and send data to Teton as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_raw" and "data_working" folders have been ignored
# using git.ignore, so links to the raw data sets are provided below.
# If you are accessing the code via GitHub, these will need to be downloaded 
# and added to a folder of the appropriate name prior to running the code.

############################
## Setup
############################

## Load packages
lapply(c("tidyverse", "cowplot","lubridate",
         "data.table","patchwork", "here"), require, character.only=T)

## Load datasets
## Import site data from Appling et al.
# https://www.sciencebase.gov/catalog/item/59bff64be4b091459a5e098b
site <- fread("data_raw/site_data.tsv")

## Import StreamLight from Savoy
# https://www.sciencebase.gov/catalog/item/5f974adfd34e198cb77db168
SL <- read.table("data_raw/StreamLight_site_information_and_parameters.txt", header=T)
colnames(SL)[colnames(SL) == "Site_ID"] <- "site_name"

## Secondary stream order source from hypoxia data set
# https://www.sciencebase.gov/catalog/item/606f60afd34ef99870188ee5
## and subsetted to Appling (PC)
hyp <- fread("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
hyp <- hyp[which(hyp$DB_Source == "PC"), c("SiteID","ORD_STRA","NHD_STREAMORDE")]
colnames(hyp)[which(colnames(hyp) == "SiteID")] <- "site_name"

## Merge StreamLight with site data, and then hypoxia dataset
df <- left_join(site, SL, "site_name")
df <- left_join(df, hyp, "site_name")

## Subset for site information regarding flow impediments, order, etc.
sub <- df[,c("site_name","long_name","StreamOrde",
             "site_type","struct.canal_flag","struct.dam_flag",
             "struct.npdes_flag","ORD_STRA","NHD_STREAMORDE")]

##################################################
## Select sites based on data quality
##################################################

#### FIRST FILTERS ####

# Import and subset model diagnostics from Appling et al.
# https://www.sciencebase.gov/catalog/item/59eb9bafe4b0026a55ffe382 
diagnostics <- read.table("data_raw/diagnostics.tsv",sep = "\t", header=T)
diagnostics <- diagnostics[which(diagnostics$site %in% sub$site_name),]

# Filtering for Rhat (Gelman-Rubin statistic) below 1.05, because this provides
# an indication for how well the streamMetabolizer output performed

# 356 sites at the start

# quick histograms to visualize the thresholds at which I'm filtering
hist(diagnostics$err_obs_iid_sigma_Rhat)
hist(diagnostics$err_proc_iid_sigma_Rhat)
hist(diagnostics$neg_GPP)
hist(diagnostics$pos_ER)

highq_sites <- diagnostics[which(
  # Instead of filtering by K600 Rhat values for the whole site, going to filter for daily
  # records using the "NWIS" dataframe below
  
  # the following two filters remove sites at which the error terms poorly converged
                                   diagnostics$err_obs_iid_sigma_Rhat < 1.05 & # 66 sites drop off
                                   diagnostics$err_proc_iid_sigma_Rhat < 1.05 & # 2 sites drops off
                                     
  # the following two filters refer to "percent of GPP estimates < -0.5 
  # and the percent of ER estimates > 0.5 (both are biologically unrealistic outcomes)"
                                   
  # at neg_GPP < 5/25/50, 48/5/2 more sites disappear

                                   diagnostics$neg_GPP < 25 &
                                     
  # at pos_ER < 5/25/50, 62/20/5 more sites disappear

                                   diagnostics$pos_ER < 25),]

# 294 records remaining after diagnostic filtering

highq_site_names <- unique(highq_sites$site) ## 261 sites at the end

# Subset s based on high quality sites and site type and flags
s <- sub[which(sub$site_name %in% highq_site_names),]

# Removing all except stream sites. Possible categories include:
# ST (stream), ST-CA (canal), ST-DCH (ditch), ST-TS (tidal stream), or SP (spring)
s <- s[which(s$site_type == "ST"),] ## 8 sites drop off

# NOTE: THERE IS NO DAM FILTER.

# Additional information from the site_data.xml file:
# "4. The fields struct.canal_flag, struct.dam_flag, struct.npdes_flag are flags for likely site suitability for metabolism modeling, based on site proximity to infrastructure that could affect metabolism estimates.
# 4a. The flagging focused on three feature types:  permitted National Pollution Discharge Elimination System point sources (see USEPA_NPDES in Source Citations), dams identified in the National Inventory of Dams (USACE_NID), and canals and ditches identified in the National Hydrography Database (NHDPlusV2). National data layers for these features of interest were clipped to watersheds upstream of the sites in our dataset. The geodetic distance between each stream monitoring site and the nearest upstream feature of each type was calculated using the function GenerateNearTable in the arcpy library in Python 3.6.
# 4b. The mean 80% oxygen turnover distance on each day (DO.tdist80) was pulled from "Daily metabolism estimates and predictors", and the 50th, 80th, and 95th percentiles of the daily 80% turnover distances were then computed for each site.
# 4c. Distances between the site and upstream features were compared with the percentiles of 80% O2 turnover distance. Site flags are numeric to indicate whether the nearest feature of each type was farther than the 95th, 80th, 50th, or 0th percentile; those four cases are represented by site flag values of 95, 80, 50, or 0, respectively, such that a value of 95 indicates the least probable interference from a structure of a given type."

# Import time series of streamMetabolizer-generated predictions
# https://www.sciencebase.gov/catalog/item/59eb9c0ae4b0026a55ffe389
NWIS <- read.table("data_raw/daily_predictions.tsv", sep='\t', header = TRUE)
NWIS$date <- ymd(NWIS$date)

# Due to date issues below, forcing the dataframe to present dates in ascending order
NWIS <- NWIS %>%
  group_by(site_name) %>%
  arrange(date) %>% # ascending is the default
  ungroup()

# Instead of filtering by site-level diagnostics above, now filtering by daily Rhat values
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

## Subset to sites in highq_sites
NWIS_sub <- NWIS_sub[which(NWIS_sub$site_name %in% s$site_name),]
# So, with diagnostics/site-specific filters revised, removed 28% of data records

# Confirm remaining number of sites
length(levels(as.factor(NWIS_sub$site_name))) ## 253 sites

#### SECOND FILTERS ####

# time-series evaluation for gaps and quantity of days

## Identify which sites have the most continuous data
NWIS_sub$doy <- yday(NWIS_sub$date)
NWIS_sub$year <- year(NWIS_sub$date)

## count days per year
dat_per_year <- NWIS_sub %>%
  count(site_name, year)

hist(dat_per_year$n)

## identify the max day gap per year
gap_per_year <- NWIS_sub %>%
  group_by(site_name, year) %>%
  mutate(gap = doy - lag(doy, default=doy[1], order_by = doy)) %>%
  ungroup()
# The code above really struggles between the dplyr/plyr packages it seems.
# So, ALWAYS double check the output above to be sure there aren't any gaps
# that are negative, as that indicates the data isn't being grouped properly.

maxgap <- gap_per_year %>%
  group_by(site_name, year) %>%
  summarize_at(.vars = "gap", .funs = max)

# Making some changes, because the shortened time frames (3 months)
# need to be included prior to filtering out for long gaps

## merge with number of days per year
sub_by_gap <- merge(maxgap, dat_per_year, by=c("site_name","year")) # 1619 site-years

# Changing to allow for season-long data to try and up the number of sites
## at least 90 days per year (minimum possible 3 month period)
sub_by_gap1 <- sub_by_gap[which(sub_by_gap$n > 90),] # 1327 site-years
sub_by_gap_sum <- sub_by_gap1 %>% group_by(site_name) %>% count() # 242 unique sites

## subset for sites with a max gap of 14 days
sub_by_gap2 <- sub_by_gap1[which(sub_by_gap1$gap <= 14),] # 784 site-years
length(levels(as.factor(sub_by_gap2$site_name))) # 207 unique sites

# Renaming dataset with sites, number of days, and gap records by year
high_q <- sub_by_gap2

## Subset NWIS_sub sites
TS <- NWIS_sub[which(NWIS_sub$site_name %in% high_q$site_name),]

## Subset to YEARS that meet criteria
sub_by_gap2$site_year <- paste(sub_by_gap2$site_name,sub_by_gap2$year,sep = "_")
TS$site_year <- paste(TS$site_name, TS$year,sep = "_")
TS <- TS[which(TS$site_year %in% sub_by_gap2$site_year),] # Site data
TS_site <- s[which(s$site_name %in% high_q$site_name),] # Site info

## Attach the median GPP
TS$GPP_temp <- TS$GPP
# Where GPP is negative, replace values with a small, randomly generated number
TS[which(TS$GPP < 0),]$GPP_temp <- sample(exp(-6):exp(-4), 1)
# Calculate site-level mean and max GPP
TS_gpp <- TS %>%
  group_by(site_name) %>%
  summarise_at(.vars = "GPP_temp", .funs = c(mean, max))
colnames(TS_gpp) <- c("site_name","GPP_mean","GPP_max")
# And append these to the larger dataset
TS_site <- left_join(TS_site, TS_gpp, by="site_name")

###########################################################################
## Choose river-years
###########################################################################

# There are 207 unique sites in the TS_site dataset, with these new filters imposed.

# first, I'm joining the full names to the TS dataset to better ID them when plotting
name_bridge <- TS_site %>%
  select(site_name, long_name)

TS <- TS %>%
  left_join(name_bridge)

# subset necessary data
site_subset <- TS # GPP data

TS_site_subset <- df[which(df$site_name %in% site_subset$site_name),] # site info

site_subset_numdays <- rbind(sub_by_gap2[which(sub_by_gap2$site_name %in% 
                                                site_subset$site_name),]) # data gaps by site-year

colnames(site_subset_numdays) <- c("site_name","year","max_gap","Ndays","site_year")

###################################################
## Check other covariate data quality
###################################################

site_sub_list <- split(site_subset, site_subset$site_name)

plotting_covar <- function(x) {
  
  df <- x
  p <- plot_grid(

    ggplot(df, aes(date, GPP))+
      geom_point(color="chartreuse4", size=2)+
      geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="darkolivegreen4")+
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=df$long_name[1])+
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

# only plotting one as a test here, but this will be handy once results come back
plotting_covar(site_sub_list$nwis_14206950)

lapply(site_sub_list, function(x) ggsave(plot = plotting_covar(x),
                                         filename = paste("figures/teton_summer22/site_covariate_plots/",x$site_name[1],"covar.jpg",sep = ""), 
                                         width = 10, 
                                         height = 8))

###########################
## Export
###########################

## NWIS site subset
saveRDS(site_subset, "data_working/NWIS_207sites_subset2022.rds") # GPP data
saveRDS(TS_site_subset, "data_working/NWIS_207sitesinfo_subset2022.rds") # Site data
saveRDS(site_subset_numdays,"data_working/NWIS_207sites_Ndays2022.rds") # Gaps in data

# End of script.
