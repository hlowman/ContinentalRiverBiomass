## Experimenting with Gaps: STAN Model Fit Script
## October 13, 2021
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to practice fitting the Ricker model to 3 sites of data
# but playing around with data missingness to examine model
# performance.

# I've commented out those steps that I feel I don't
# need to perform. I've also changed the filepaths to match my
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

# Import and subset model diagnostics
# https://www.sciencebase.gov/catalog/item/59eb9bafe4b0026a55ffe382
diagnostics <- read.table("data_raw/diagnostics.tsv",sep = "\t", header=T)
diagnostics <- diagnostics[which(diagnostics$site %in% sub$site_name),]
highq_sites <- diagnostics[which(diagnostics$K600_daily_sigma_Rhat < 1.05 & 
                                   diagnostics$err_obs_iid_sigma_Rhat < 1.05 &
                                   diagnostics$err_proc_iid_sigma_Rhat < 1.05 &
                                   diagnostics$neg_GPP < 15 &
                                   diagnostics$pos_ER < 15),] #229
highq_site_names <- unique(highq_sites$site) ## 208

# Subset s based on high sites and site type and flags
s <- sub[which(sub$site_name %in% highq_site_names),] ## 208
s <- s[which(s$site_type == "ST"),] ## 204
s <- s[which(s$struct.dam_flag %in% c(NA,"95")),] ## 82
# which have light from Phil
s_l <- s[!is.na(s$StreamOrde),] # 36

# Import time series
NWIS <- read.table("data_raw/daily_predictions.tsv", sep='\t', header = TRUE)
NWIS$date <- as.POSIXct(as.character(NWIS$date), format="%Y-%m-%d")

## Subset columns and sites
NWIS_sub <- NWIS[,c("site_name","date","GPP","GPP.lower","GPP.upper", 
                    "GPP.Rhat","ER","ER.lower","ER.upper","K600",
                    "K600.lower","K600.upper","temp.water",
                    "discharge","shortwave","velocity")]
colnames(NWIS_sub) <- c("site_name","date","GPP","GPP.lower","GPP.upper", 
                        "GPP.Rhat","ER","ER.lower","ER.upper","K600",
                        "K600.lower","K600.upper","temp","Q","light",
                        "velocity")

## Subset to sites in high_sites (sites with high confidence rating and limited dam interference)
NWIS_sub <- NWIS_sub[which(NWIS_sub$site_name %in% s$site_name),]
# Confirm
length(levels(as.factor(NWIS_sub$site_name))) ## 82 when subsetting for s

## Identify which sites have the most continuous data
NWIS_sub$doy <- yday(NWIS_sub$date)
NWIS_sub$year <- year(NWIS_sub$date)

## count days per year
dat_per_year <- NWIS_sub %>%
  group_by(site_name, year) %>%
  count()

## identify the max day gap per year
gap_per_year <- NWIS_sub %>%
  group_by(site_name, year) %>%
  mutate(gap = doy - lag(doy, default=doy[1]))

maxgap <- gap_per_year %>%
  group_by(site_name, year) %>%
  summarize_at(.vars = "gap", .funs = max)

## subset for sites with a max gap of 14 days
sub_by_gap <- maxgap[which(maxgap$gap <= 14),]
length(levels(as.factor(sub_by_gap$site_name))) #77

## merge with number of days per year
sub_by_gap <- merge(sub_by_gap, dat_per_year, by=c("site_name","year"))

## at least 275 days per year
sub_by_gap <- sub_by_gap[which(sub_by_gap$n > 275),]
sub_by_gap_sum <- sub_by_gap %>% group_by(site_name) %>% count()
high_q <- sub_by_gap_sum[which(sub_by_gap_sum$n >= 2),]

## Subset NWIS_sub
TS <- NWIS_sub[which(NWIS_sub$site_name %in% high_q$site_name),]

## Subset to years that meet criteria
sub_by_gap$site_year <- paste(sub_by_gap$site_name,sub_by_gap$year,sep = "_")
TS$site_year <- paste(TS$site_name, TS$year,sep = "_")
TS <- TS[which(TS$site_year %in% sub_by_gap$site_year),]
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
## Choose 3 sites
###########################################################################

# There are 34 unique sites in the TS_site dataset, so I am going to use these.

# first, I'm joining the full names to the TS dataset to better ID them when plotting
name_bridge <- TS_site %>%
  select(site_name, long_name)

TS <- TS %>%
  left_join(name_bridge)

# count all site-years of GPP data
site_yr_count <- TS %>%
  mutate(yearf = factor(year)) %>%
  group_by(site_name, yearf) %>%
  summarize(meanGPP = mean(GPP)) %>%
  ungroup() %>%
  count(site_name) %>%
  mutate(n_f = factor(n))

# chose sites
chosen_sites <- c("nwis_02266300", "nwis_14206950", "nwis_01608500")

# filter dataset only for selected sites
TS3 <- TS %>%
  filter(site_name %in% chosen_sites)

# subset necessary data
site_subset <- TS3 # GPP data

TS_site_subset <- df[which(df$site_name %in% site_subset$site_name),] # site info

site_subset_numdays <- rbind(sub_by_gap[which(sub_by_gap$site_name %in% 
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

# plotting them all
plotting_covar(site_sub_list$nwis_02266300) # remove 3yrs
plotting_covar(site_sub_list$nwis_14206950) # remove 2yrs
plotting_covar(site_sub_list$nwis_01608500) # remove 1yr

###########################
## Export
###########################

## NWIS site subset
saveRDS(site_subset, "data_working/NWIS_3sites_subset.rds") # GPP data
saveRDS(TS_site_subset, "data_working/NWIS_3sitesinfo_subset.rds") # Site data
saveRDS(site_subset_numdays,"data_working/NWIS_3sites_Ndays.rds") # Gaps in data

# End of script.
