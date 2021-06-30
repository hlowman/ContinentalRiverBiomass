## More Practice STAN Model Fit Script
## June 23, 2021
## Heili Lowman

# I'll be modifying some of Joanna's code from the RiverBiomass repository
# to practice fitting the Ricker model to 1 year of data at a "good" and
# "bad" site.

# I've commented out those steps that I feel, for the time being, I don't
# need to perform, and I've changed the appropriate filepaths to match my
# repository structure.

## Subset data from hypoxia database that is already linked to NHD
## JRB

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork", "here"), require, character.only=T)

############################
## To create linked file
############################

## Import site data from Appling
# setwd("../data")
site <- fread("data_raw/site_data.tsv")
names(site)

## Import StreamLight from Savoy
# https://www.sciencebase.gov/catalog/item/5f974adfd34e198cb77db168
SL <- read.table("data_raw/StreamLight_site_information_and_parameters.txt", header=T)
colnames(SL)[colnames(SL) == "Site_ID"] <- "site_name"

## Secondary stream order source from hypoxia data set
#https://www.sciencebase.gov/catalog/item/606f60afd34ef99870188ee5
## and subsetted to Appling (PC)
hyp <- fread("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
hyp <- hyp[which(hyp$DB_Source == "PC"), c("SiteID","ORD_STRA","NHD_STREAMORDE")]
colnames(hyp)[which(colnames(hyp) == "SiteID")] <- "site_name"

## Merge
df <- left_join(site, SL, "site_name")
df <- left_join(df, hyp, "site_name")

colnames(df)

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
head(NWIS)

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
length(levels(as.factor(sub_by_gap$site_name)))
## merge with number of days per year
sub_by_gap <- merge(sub_by_gap, dat_per_year, by=c("site_name","year"))
## at least 275 days per year
sub_by_gap <- sub_by_gap[which(sub_by_gap$n > 275),]
sub_by_gap_sum <- sub_by_gap %>% group_by(site_name) %>% count()
high_q <- sub_by_gap_sum[which(sub_by_gap_sum$n >= 2),]

## Subset NWIS_sub
TS <- NWIS_sub[which(NWIS_sub$site_name %in% high_q$site_name),]
## only sites with two or more years

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
## Choose one river year from a "good" river
###########################################################################

# looking for a medium-sized stream with appropriate-seeming
# responses in GPP to disturbance
## plot
# SOUTH BRANCH POTOMAC RIVER NEAR SPRINGFIELD, WV
sid <- "nwis_01608500"
years <- c(2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018)
TS_site[which(TS_site$site_name == sid),]

plot_grid(
  ggplot(TS[which(TS$site_name == sid),], aes(date, GPP_temp))+
    geom_line()+labs(title=TS_site[which(TS_site$site_name == sid),]$long_name),
  ggplot(TS[which(TS$site_name == sid & TS$year %in% years),], aes(date, GPP_temp))+
    geom_line(),
  ncol = 1)

## plot my chosen site
ggplot(TS[which(TS$site_name == "nwis_01608500"),], 
       aes(date, GPP))+
  geom_line() +
  labs(title = "South Branch Potomac River, Springfield, WV")
# something again appears to be happening to make values negative
# so going to use 2015 going forward

## "good": nwis_01608500 2015 (South Branch Potomac, WV)

site_subset <- rbind(TS[which(TS$site_name == "nwis_01608500" & 
                                TS$year %in% c(2015)),])

TS_site_subset <- df[which(df$site_name %in% site_subset$site_name),]

## Save sub_by_gap info
site_subset_numdays <- rbind(sub_by_gap[which(sub_by_gap$site_name == "nwis_01608500" & sub_by_gap$year %in% c(2015)),])
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

plotting_covar(site_sub_list$nwis_01608500)

lapply(site_sub_list, function(x) ggsave(plot = plotting_covar(x),filename = paste("figures/site_covariate_plots/",x$site_name[1],"covar.jpg",sep = ""), width = 8, height = 6))


###########################
## Export
###########################

## NWIS site subset
saveRDS(site_subset, "data_working/NWIS_1site_subset_good_2015.rds")
saveRDS(TS_site_subset, "data_working/NWIS_1siteinfo_subset_good_2015.rds")
saveRDS(site_subset_numdays,"data_working/NWIS_1site_Ndays_good_2015.rds")

#### Stopped here during test run on June 23, 2021.

##############################
## Plot for talk
###########################

# df <- site_sub_list$nwis_14211010
# ratio_QL <- max(df$light)/max(df$Q)
# GPP_plot <- ggplot(df, aes(date, GPP))+
#   geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="chartreuse4")+
#       geom_point(color="chartreuse4", size=2)+geom_line(color="chartreuse4", size=1)+
#       labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'))+
#       theme(legend.position = "none",
#             panel.background = element_rect(color = "black", fill=NA, size=1),
#             axis.title.x = element_blank(), axis.text.x = element_blank(),
#             axis.text.y = element_text(size=12),
#             axis.title.y = element_text(size=12))+
#   geom_vline(xintercept = as.POSIXct("2013-01-01", format="%Y-%m-%d"),size=1, linetype="dashed")
     
# data_plot <- ggplot(df, aes(date, Q*ratio_QL))+
#       geom_point(data=df, aes(date, light), size=1.5, color="darkgoldenrod3")+
#       geom_line(size=1, color="deepskyblue4")+
#       scale_y_continuous(sec.axis = sec_axis(~./ratio_QL, name=expression("Daily Q (cms)")))+
#       labs(y=expression('Daily PPFD'))+# ('*~mu~mol~ m^-2~d^-1*')'), x="Date")+
#       theme(legend.position = "none",
#             panel.background = element_rect(color = "black", fill=NA, size=1),
#             axis.title.x = element_blank(), axis.text = element_text(size=12),
#             axis.title.y.left = element_text(size=12, color="darkgoldenrod3"),
#             axis.title.y.right = element_text(size=12, color="deepskyblue4"),
#             axis.text.x = element_text(angle=25, hjust = 1),
#             strip.background = element_rect(fill="white", color="black"),
#             strip.text = element_text(size=15))+
#   geom_vline(xintercept = as.POSIXct("2013-01-01", format="%Y-%m-%d"),size=1, linetype="dashed")
     
# GPP_plot + data_plot + plot_layout(ncol = 1)
 
# ggplot(df, aes(light, GPP))+
#   geom_point(size=1.5)+
#   labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x="Daily PPFD")+
#   theme_bw()+
#   theme(legend.position = "none",
#         panel.background = element_rect(color = "black", fill=NA, size=1),
#         axis.text = element_text(size=14),
#         axis.title = element_text(size=18))

# End of script.
