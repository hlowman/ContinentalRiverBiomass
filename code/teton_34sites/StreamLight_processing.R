## StreamLight Output Script
## July 13, 2021
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to practice fitting the Ricker model to ~30 site-years of data.
# The compiled data will then be run by sending the job to Teton (UWyoming).

# I've commented out those steps that I feel, for the time being, I don't
# need to perform, and I've changed the appropriate filepaths to match my
# repository structure.

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork", "here"), require, character.only=T)

# Load in datasets created in NWIS_RiverSelection script.
NWIS_34sites_subset <- readRDS("data_working/NWIS_34sites_subset.rds")
NWIS_34sitesinfo_subset <- readRDS("data_working/NWIS_34sitesinfo_subset.rds")
NWIS_34sites_Ndays <- readRDS("data_working/NWIS_34sites_Ndays.rds")

## Site of interest that has StreamLight data
sites <- NWIS_34sitesinfo_subset$site_name
sites_files <- rep(NA, length(sites))
for(i in 1:length(sites)){
  sites_files[i] <- paste(sites[i],"_input_output.txt", sep="")
}

# need to identify the intersection of files with available light data
filenames <- intersect(sites_files, list.files(here("data_raw", "individual_files")))

## Import StreamLight
SL <- ldply(filenames, function(filename) {
  d <- read.table(here("data_raw", "individual_files", filename), header = T, sep = "\t")
  d$file <- filename
  return(d)
})

## take the mean daily incoming PAR at the surface
SL_split <- split(SL, SL$file)

meandaily_PAR <- function(y){
  df <- y %>%
  group_by(jday) %>%
  summarize_at(.vars = c("DOY","Year","PAR_surface","PAR_turb"), .funs = mean, na.rm = TRUE)
  
  df$origin <- as.Date(paste0(df$Year, "-01-01"),tz = "UTC") - days(1)
  df$Date <- as.Date(df$DOY, origin = df$origin, tz = "UTC") 
  
  return(df)
}

SL_daily <- lapply(SL_split, function(x) meandaily_PAR(x))

## visualize
ggplot(SL_daily$nwis_14206950_input_output.txt, aes(Date, PAR_surface))+
  geom_point()+
  geom_point(aes(Date, PAR_turb), color="blue") +
  labs(x = "Date",
       y = "Daily PAR",
       title = "Fanno Creek - Durham, OR")

###########################
## Subset and evaluate NA
###########################
SF_df <- ldply(SL_daily, data.frame)

## add site_name
SF_df$site_name <- substr(SF_df$.id, 1, nchar(SF_df$.id)-17)

site_subset <- SF_df

site_subset_split <- split(site_subset, site_subset$.id)

# what are these doing??
lapply(site_subset_split, function(x) sum(is.na(x$PAR_surface)))
lapply(site_subset_split, function(x) sum(is.na(x$PAR_turb)))

## Save

saveRDS(site_subset, "data_working/StreamLight_daily_34sites.rds")

# End of script.
