## StreamLight Output Script
## June 23, 2021
## Heili Lowman

# I'll be modifying some of Joanna's code from the RiverBiomass repository
# to practice fitting the Ricker model to 1 year of data at a "good" and
# "bad" site.

# I've commented out those steps that I feel, for the time being, I don't
# need to perform, and I've changed the appropriate filepaths to match my
# repository structure.

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork", "here"), require, character.only=T)

## Site of interest that has StreamLight data
sites <- c("nwis_01608500")
sites_files <- rep(NA, length(sites))
for(i in 1:length(sites)){
  sites_files[i] <- paste(sites[i],"_input_output.txt", sep="")
}

## Import streamlight

SL <- ldply(sites_files, function(filename) {
  d <- read.table(here("data_raw", "individual_files", filename), header = T, sep = "\t")
  d$file <- filename
  return(d)
})

## take the mean daily incoming PAR at the surface
SL_split <- split(SL, SL$file)
View(SL_split$nwis_01608500_input_output.txt)

meandaily_PAR <- function(y){
  df <- y %>%
  group_by(jday) %>%
  summarize_at(.vars = c("DOY","Year","PAR_surface","PAR_turb"), .funs = mean, na.rm = TRUE)
  
  df$origin <- as.Date(paste0(df$Year, "-01-01"),tz = "UTC") - days(1)
  df$Date <- as.Date(df$DOY, origin = df$origin, tz = "UTC") 
  
  return(df)
}

SL_daily <- lapply(SL_split, function(x) meandaily_PAR(x))
tail(SL_daily$nwis_01608500_input_output.txt)

## visualize
ggplot(SL_daily$nwis_01608500_input_output.txt, aes(Date, PAR_surface))+
  geom_point()+
  geom_point(aes(Date, PAR_turb), color="blue")+
  labs(title = "South Branch Potomac River, WV")


###########################
## Subset and evaluate NA
###########################
SF_df <- ldply(SL_daily, data.frame)

## add site_name
SF_df$site_name <- substr(SF_df$.id, 1, nchar(SF_df$.id)-17)

## From NWIS_RiverSelection
site_subset <- rbind(SF_df[which(SF_df$site_name == "nwis_01608500" & SF_df$Year %in% c(2015)),])

site_subset_split <- split(site_subset, site_subset$.id)

lapply(site_subset_split, function(x) sum(is.na(x$PAR_surface))) # all 0
lapply(site_subset_split, function(x) sum(is.na(x$PAR_turb))) # 60

## Save

saveRDS(site_subset, "data_working/StreamLight_daily_1riv_good_2015.rds")

# End of script.
