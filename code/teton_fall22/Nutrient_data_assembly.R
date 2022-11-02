## Nutrients data source
## November 2, 2022
## Heili Lowman

# This script will download available nutrient data from the USGS using their
# dataRetrieval package.

#### Setup ####

# Load packages
lapply(c("tidyverse", "lubridate", "here",
         "dataRetrieval"), require, character.only=T)

# Load in original list of data fed into the model.
dat_in <- readRDS("data_working/list_182sites_Qmaxnorm_allSL.rds")

#### Parameter Code Aggregation ####

# First, need to identify the parameter codes for my nutrients of interest
# I know that broadly I'll need phosphorus, nitrogen, nitrate, and ammonia,
# as listed in L. Koenig's data pull workflow.

# Load in all available codes.
pcode <- readNWISpCode("all")

# Pull out P codes
phosCds <- pcode[grep("phosphorus",
                      pcode$parameter_nm,
                      ignore.case=TRUE),]

# Pull out N codes
nitroCds <- pcode[grep("nitrogen",
                      pcode$parameter_nm,
                      ignore.case=TRUE),]

# Pull out NO3 codes
nitraCds <- pcode[grep("nitrate",
                       pcode$parameter_nm,
                       ignore.case=TRUE),]

# Pull out NH3 codes
ammoCds <- pcode[grep("ammonia",
                       pcode$parameter_nm,
                       ignore.case=TRUE),]

# Join these together
nutCds1 <- full_join(phosCds, nitroCds)
nutCds2 <- full_join(nutCds1, nitraCds)
nutCds <- full_join(nutCds2, ammoCds)

#### Test Workflow ####

# Use the NEW water quality function to test pull a single type of data 
# at a single site.
# More info here: https://cran.r-project.org/web/packages/dataRetrieval/vignettes/qwdata_changes.html
# Info about deprecated functions here: https://waterdata.usgs.gov/blog/dataretrieval/
siteNo <- "02203900"
pCode <- "00662"
start.date <- "2007-01-01"
end.date <- "2017-12-31"

wqpTest <- readWQPqw(paste0("USGS-", siteNo), # new nomenclature for sites
                     pCode,
                     startDate = start.date,
                     endDate = end.date)

# And testing for all parameter codes.
all_codes <- nutCds$parameter_cd

wqpTest2 <- readWQPqw(paste0("USGS-", siteNo), 
                     all_codes,
                     startDate = start.date,
                     endDate = end.date)

#### USGS/WQP Data Pull ####

# Create new list of sites for which we need to pull data.
my_sites <- names(dat_in)
my_sites <- str_replace_all(my_sites, 'nwis_', '') # removes leading characters
my_sites <- paste0("USGS-", my_sites) # adds correct leading characters

# Large data pull.
wqpNuts <- readWQPqw(my_sites, 
                     all_codes,
                     startDate = start.date,
                     endDate = end.date)

# And export raw data pull.
saveRDS(wqpNuts, "data_raw/USGS_WQP_Nutrients_2022_11_02.rds")

#### Summarizing data ####

# For now, summarizing data by site since data coverage may be all over the
# place.

summNuts <- wqpNuts %>%
  group_by(MonitoringLocationIdentifier, 
           CharacteristicName,
           ResultMeasure.MeasureUnitCode) %>%
  summarize(mean_value = mean(ResultMeasureValue, na.rm=TRUE)) %>%
  ungroup()

# Examine data
no3 <- summNuts %>%
  filter(ResultMeasure.MeasureUnitCode == "mg/l asNO3")

hist(no3$mean_value) # concentrations generally seem pretty low.

p <- summNuts %>%
  filter(CharacteristicName == "Phosphorus")

hist(p$mean_value) # same.

# End of script.
