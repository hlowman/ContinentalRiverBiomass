## Nutrients data source
## November 2, 2022
## Heili Lowman

# This script will download available nutrient data from the USGS using their
# dataRetrieval package.

#### Setup ####

# Load packages
lapply(c("tidyverse", "lubridate", "here",
         "dataRetrieval", "viridis"), require, character.only=T)

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
  summarize(mean_value = mean(ResultMeasureValue, na.rm=TRUE),
            measure_count = sum(!is.na(ResultMeasureValue))) %>% # bc otherwise
  # n() function counts NAs, which we don't want
  ungroup()

# Examine data
no3 <- summNuts %>%
  filter(ResultMeasure.MeasureUnitCode == "mg/l asNO3")

hist(no3$mean_value) # concentrations generally seem pretty low.

p <- summNuts %>%
  filter(CharacteristicName == "Phosphorus")

hist(p$mean_value) # same.

#### Quick Figures ####

# Need to make columns of the analytes.
# But first, need to select analytes/columns of interest
summ_No3_filter <- summNuts %>%
  filter(ResultMeasure.MeasureUnitCode == "mg/l asNO3") %>%
  select(MonitoringLocationIdentifier, mean_value) %>%
  rename("NO3_mg_L" = "mean_value")

summ_oP_filter <- summNuts %>%
  filter(ResultMeasure.MeasureUnitCode == "mg/l as P") %>%
  filter(CharacteristicName == "Orthophosphate") %>%
  select(MonitoringLocationIdentifier, mean_value) %>%
  rename("oP_mg_L_P" = "mean_value")

# And then join back together.
summ_No3_oP <- full_join(summ_No3_filter, summ_oP_filter, by = c("MonitoringLocationIdentifier"))

(fig1 <- ggplot(summ_No3_oP, aes(x = log10(NO3_mg_L), y = log10(oP_mg_L_P))) +
  geom_point(size = 3, alpha = 0.8, aes(color = NO3_mg_L)) +
  scale_color_viridis() +
  labs(x = "Log of Nitrate (mg/L as NO3)",
       y = "Log of orthoPhosphate (mg/L as P)",
       color = "[NO3] (mg/L)",
       title = "USGS Water Quality Data (n = 114 sites)") +
  theme_bw())

# Export figure.
# ggsave(fig1,
#        filename = "figures/teton_fall22/nuts_exploration_fig.jpg",
#        width = 15,
#        height = 10,
#        units = "cm")

# Another figure to look at frequency of certain measurements by site.
(fig2 <- ggplot(summNuts, aes(x = CharacteristicName, y = measure_count)) +
    geom_boxplot(aes(color = CharacteristicName)) +
    scale_color_viridis(discrete = TRUE) +
    labs(x = "Analytes",
         y = "Count of Measurements") +
    theme_bw() +
    theme(legend.position = "none") +
    coord_flip())

# Export figure.
# ggsave(fig2,
#        filename = "figures/teton_fall22/nuts_availability_fig.jpg",
#        width = 20,
#        height = 10,
#        units = "cm")

# Looks like phosphorus and nitrate may be most abundant, but Joanna suggested
# using orthophosphate so going to move ahead with those two measures.

#### Orthophosphate ####

# Need to examine what units are present for P.
orthoP <- wqpNuts %>%
  filter(CharacteristicName == "Orthophosphate")

unique(orthoP$ResultMeasure.MeasureUnitCode)  # only one code - "mg/l as P"

# Note for future me - "Phosphorus" has ~4 codes, so that would make conversion
# a little more tricky.

# Just for another quick gut check
length(unique(orthoP$MonitoringLocationIdentifier)) # 112 sites with oP data

# Note, will need to remove measurements with NA as unit.

# FINAL UNITS: mg/L as PO4-P

#### Nitrate ####
nitrate <- wqpNuts %>%
  filter(CharacteristicName == "Nitrate")
  
unique(nitrate$ResultMeasure.MeasureUnitCode) # two codes - "mg/l asNO3", 
# and "mg/l as N", so will need to do some conversions
  
# Also need to remove measurements with NA as unit.

# Just for another quick gut check
length(unique(nitrate$MonitoringLocationIdentifier)) # 105 sites with NO3 data

# Need to convert mg/L as N to mg/L as NO3.
# 1 mg/L No3-N = 4.43 mg/L NO3
# 0.226 mg/L NO3-N = 1 mg/L NO3
wqpNuts <- wqpNuts %>%
  mutate(ResultMeasureValue_conv = case_when(CharacteristicName == "Nitrate" &
                                               ResultMeasure.MeasureUnitCode == "mg/l asNO3" ~ ResultMeasureValue*0.226,
                                             TRUE ~ ResultMeasureValue))

# FINAL UNITS: mg/L as NO3-N

#### Summary by site ####

# Summarize Orthophosphate and Nitrate data by site.
wqp_oP_NO3 <- wqpNuts %>%
  filter(CharacteristicName %in% c("Orthophosphate", "Nitrate")) %>%
  group_by(MonitoringLocationIdentifier, CharacteristicName) %>%
  summarize(mean_mg_L = mean(ResultMeasureValue_conv, na.rm = TRUE)) %>%
  ungroup()

# And export for use in other script.
saveRDS(wqp_oP_NO3, "data_working/USGS_WQP_nuts_aggsite_110222.rds")
  
# End of script.
