## Step ONE in Metabolism Modeling Workflow
## November 16, 2021
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to fit the Ricker model to a series of selected sites.

# Specifically, the workflow in this folder will work on pooling data,
# working across gaps in time series as well as models w/ and w/o the P term.

# I've commented out those steps that I feel, for the time being, I don't
# need to perform, and I've changed the appropriate filepaths to match my
# repository structure.

# Additional note: the "data_raw" and "data_working" folders have been ignored
# in all iterations of the code (on various computers), so links to the raw
# datasets are provided below and will need to be downloaded and added to a folder
# of the appropriate name prior to running the code.

## Subset data from hypoxia database that is already linked to NHD
## JRB

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork", "here"), require, character.only=T)

############################
## To create linked file
############################

## Load datasets
## Import site data from Appling
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

## Merge
df <- left_join(site, SL, "site_name")
df <- left_join(df, hyp, "site_name")

## subset
sub <- df[,c("site_name","long_name","StreamOrde",
             "site_type","struct.canal_flag","struct.dam_flag",
             "struct.npdes_flag","ORD_STRA","NHD_STREAMORDE")]

##################################################
## Select sites based on data quality
##################################################

# Using code developed in November 2021, so the newer filtering regime is imposed.

#### FIRST FILTERS ####
# diagnostics from Appling et al.

# Import and subset model diagnostics
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

highq_sites <- diagnostics[which(#diagnostics$K600_daily_sigma_Rhat < 1.05 & 
  # Instead of filtering by K600 Rhat values for the whole site, going to filter for daily
  # records using the "NWIS" dataframe below
  
                                   diagnostics$err_obs_iid_sigma_Rhat < 1.05 & # 66 sites drop off
                                   diagnostics$err_proc_iid_sigma_Rhat < 1.05 & # 2 sites drops off
                                     
                                   # the following two filters refer to "percent of GPP estimates < -0.5 
                                   # and the percent of ER estimates > 0.5 
                                   # (both are biologically unrealistic outcomes)"
                                   # at neg_GPP < 5, 48 more sites disappear
                                   # at neg_GPP < 25, 5 sites disappear
                                   # at neg_GPP < 50, only 2 sites drop off
                                   diagnostics$neg_GPP < 25 &
                                     
                                   # at pos_ER < 5, 62 more sites disappear
                                   # at pos_ER < 25, 20 sites disappear
                                   # at pos_ER < 50, 5 sites drop off
                                   diagnostics$pos_ER < 25),]

# 294 records remaining after diagnostic filtering

highq_site_names <- unique(highq_sites$site) ## 294 sites at the end

# Subset s based on high sites and site type and flags
s <- sub[which(sub$site_name %in% highq_site_names),]

# Removing all except stream sites. Possible categories include:
# ST (stream), ST-CA (canal), ST-DCH (ditch), ST-TS (tidal stream), or SP (spring)
s <- s[which(s$site_type == "ST"),] ## 8 sites drop off

# NOTE: REMOVING THE DAM FILTER. This will instead be addressed with variable application of the
# "P" term in future model fits.

# Removing possible interference from dams. From the Appling site_data.xml file:
# "a value of 95 indicates the least probable interference from a structure of a given type"
# s <- s[which(s$struct.dam_flag %in% c(NA,"95")),] ## 122 sites drop off, 82 left
#s_test <-  s[which(s$struct.dam_flag %in% c(NA,"80", "95")),] ## 100 sites drop off, 104 left
#s_test2 <-  s[which(s$struct.dam_flag %in% c(NA, "50", "80", "95")),] ## 81 sites drop off, 123 left

# Additional information from the site_data.xml file:
# "4. The fields struct.canal_flag, struct.dam_flag, struct.npdes_flag are flags for likely site suitability for metabolism modeling, based on site proximity to infrastructure that could affect metabolism estimates.
# 4a. The flagging focused on three feature types:  permitted National Pollution Discharge Elimination System point sources (see USEPA_NPDES in Source Citations), dams identified in the National Inventory of Dams (USACE_NID), and canals and ditches identified in the National Hydrography Database (NHDPlusV2). National data layers for these features of interest were clipped to watersheds upstream of the sites in our dataset. The geodetic distance between each stream monitoring site and the nearest upstream feature of each type was calculated using the function GenerateNearTable in the arcpy library in Python 3.6.
# 4b. The mean 80% oxygen turnover distance on each day (DO.tdist80) was pulled from "Daily metabolism estimates and predictors", and the 50th, 80th, and 95th percentiles of the daily 80% turnover distances were then computed for each site.
# 4c. Distances between the site and upstream features were compared with the percentiles of 80% O2 turnover distance. Site flags are numeric to indicate whether the nearest feature of each type was farther than the 95th, 80th, 50th, or 0th percentile; those four cases are represented by site flag values of 95, 80, 50, or 0, respectively, such that a value of 95 indicates the least probable interference from a structure of a given type."

# which have light from Phil, this isn't so much a filter as a possibility
# s_l <- s[!is.na(s$StreamOrde),] # 36 sites have StreamLight data

# Import time series of streamMetabolizer-generated predictions
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
  filter(K600.Rhat < 1.05) # same

# Removed 5% of days (470,763 records remaining of the original 490,907)

## Subset columns and sites
NWIS_sub <- NWIS_ed[,c("site_name","date","GPP","GPP.lower","GPP.upper", 
                    "GPP.Rhat","ER","ER.lower","ER.upper","K600",
                    "K600.lower","K600.upper","temp.water",
                    "discharge","shortwave","velocity")]
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
  group_by(site_name, year) %>%
  count()

# and visualize distribution of days per year across all site-years
hist(dat_per_year$n)

## identify the max day gap per year
gap_per_year <- NWIS_sub %>%
  group_by(site_name, year) %>%
  mutate(gap = doy - lag(doy, default=doy[1], order_by = doy))
# now, no more negative values in the gap column - woohoo!

# So, the following code shows that the dates are now arranged in ascending order
# When earlier this was not the case using the same code
test <- gap_per_year %>% filter(site_name == "nwis_040871488" & year == 2011) # YAY! It's working.
# Switching the above code to lubridate's ymd() function alone didn't seem to help
# so I'm going to force it to be in ascending order above using arrange().

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

# Commenting out the below, since having more than 1 year of data doesn't matter
# for our purposes
# high_q <- sub_by_gap_sum[which(sub_by_gap_sum$n >= 2),]

high_q <- sub_by_gap2

## Subset NWIS_sub
TS <- NWIS_sub[which(NWIS_sub$site_name %in% high_q$site_name),]

## Subset to years that meet criteria
sub_by_gap2$site_year <- paste(sub_by_gap2$site_name,sub_by_gap2$year,sep = "_")
TS$site_year <- paste(TS$site_name, TS$year,sep = "_")
TS <- TS[which(TS$site_year %in% sub_by_gap2$site_year),]
TS_site <- s[which(s$site_name %in% high_q$site_name),]

## Attach the median GPP
TS$GPP_temp <- TS$GPP
TS[which(TS$GPP < 0),]$GPP_temp <- sample(exp(-6):exp(-4), 1)
TS_gpp <- TS %>%
  group_by(site_name) %>%
  summarise_at(.vars = "GPP_temp", .funs = c(mean, max))
colnames(TS_gpp) <- c("site_name","GPP_mean","GPP_max")
TS_site <- left_join(TS_site, TS_gpp, by="site_name")

## Assign a stream order classification
TS_site$order_group <- "NA"
TS_site[which(TS_site$NHD_STREAMORDE %in% c(1,2)),]$order_group <- "small"
TS_site[which(TS_site$NHD_STREAMORDE %in% c(3,4,5)),]$order_group <- "mid"
TS_site[which(TS_site$NHD_STREAMORDE >= 6),]$order_group <- "large"

###########################################################################
## Choose river-years
###########################################################################

# There are 207 unique sites in the TS_site dataset, with these new filters imposed.

# first, I'm joining the full names to the TS dataset to better ID them when plotting
name_bridge <- TS_site %>%
  select(site_name, long_name)

TS <- TS %>%
  left_join(name_bridge)

# plot all site-years of GPP data
ggplot(TS, aes(date, GPP_temp)) +
  geom_line() +
  labs(x = "Date",
       y = "GPP",
       title = "Site-Years for Second Teton Job") +
  facet_wrap(.~long_name, scales = "free")

# plot coefficient of variation (sd/mean) of discharge by site
(fig_cvQ <- TS %>%
  group_by(site_name) %>%
  summarize(cvQ = (sd(Q, na.rm = TRUE)/mean(Q, na.rm = TRUE))) %>%
  ungroup() %>%
  left_join(s, by = "site_name") %>%
  select(site_name, long_name, cvQ, struct.dam_flag) %>%
  ggplot(aes(cvQ, struct.dam_flag)) +
  geom_point() +
  theme_bw() +
  labs(x = "Coefficient of Variation of Discharge (Q)",
       y = "Probable Interference from Dams (95 = least probable)",
       title = "207 Sites in Second Teton Job"))

(fig_cvQ2 <- TS %>%
    group_by(site_name) %>%
    summarize(cvQ = (sd(Q, na.rm = TRUE)/mean(Q, na.rm = TRUE))) %>%
    ungroup() %>%
    left_join(s, by = "site_name") %>%
    select(site_name, long_name, cvQ, NHD_STREAMORDE, struct.dam_flag) %>%
    mutate(so = factor(NHD_STREAMORDE)) %>%
    ggplot(aes(cvQ, so)) +
    geom_point(aes(color = struct.dam_flag), size = 3) +
    geom_boxplot(fill = NA) +
    #scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    theme_bw() +
    labs(x = "Coefficient of Variation of Discharge (Q)",
         y = "Stream Order",
         color = "Probable Interference from\nDams (95 = least probable)",
         title = "207 Sites in Second Teton Job"))

# Larger streams (higher order) tend to have greater probability of interference from dams
# Larger streams also tend to have lower coefficients of variation

# site with multiple time segments
# plot Allegheny River data
TS %>%
  filter(long_name == "Allegheny River at Franklin, PA") %>%
  ggplot(aes(date, GPP_temp)) +
  geom_line() +
  labs(x = "Date",
       y = "GPP")

# site with poor GPP response to Q
# plot Reedy Creek data
TS %>%
  filter(long_name == "REEDY CREEK NEAR VINELAND, FL") %>%
  ggplot(aes(date, GPP_temp)) +
  geom_line() +
  labs(x = "Date",
       y = "GPP")

# site with both multiple time segments and poor GPP response to Q
# plot Santa Margarita River data
TS %>%
  filter(long_name == "SANTA MARGARITA R NR TEMECULA CA") %>%
  ggplot(aes(date, GPP_temp)) +
  geom_line() +
  labs(x = "Date",
       y = "GPP")

# subset necessary data
site_subset <- TS %>% # GPP data
  filter(site_name %in% c("nwis_03025500", "nwis_02266300", "nwis_11044000"))

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
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=df$site_name[1])+
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

lapply(site_sub_list, function(x) ggsave(plot = plotting_covar(x),filename = paste("figures/pooling_practice/site_covariate_plots/",x$site_name[1],"covar.jpg",sep = ""), width = 8, height = 6))


###########################
## Export
###########################

## NWIS site subset
saveRDS(site_subset, "data_working/NWIS_pooling_3sites_subset.rds") # GPP data
saveRDS(TS_site_subset, "data_working/NWIS_pooling_3sitesinfo_subset.rds") # Site data
saveRDS(site_subset_numdays,"data_working/NWIS_pooling_3sites_Ndays.rds") # Gaps in data

# End of script.
