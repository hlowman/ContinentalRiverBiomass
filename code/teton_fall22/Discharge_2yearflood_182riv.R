## Estimating 10 year floods
## October 17, 2022
## Heili Lowman

# This code has been adapted from Joanna's code found at:
# https://github.com/jrblaszczak/RiverBiomass/blob/main/code/Discharge_TS_2yearflood_6riv.R

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork","dataRetrieval"), require, character.only=T)

# Read in main dataframe so that we can extract site names
df <- readRDS("data_working/teton_182rivers_model_diags_101522.rds")

# Site numbers
site_numbers <- unique(df$site_name)

# And remove the "nwis_" from each
site_numbers <- str_remove_all(site_numbers, "[nwis_]")

# Site information
NWIS_Info <- lapply(site_numbers, function(x) dataRetrieval::readNWISsite(x))

# Mean daily discharge
parameterCd <- "00060"

# Raw daily data for the past 50 years
rawDailyData <- lapply(site_numbers, function(y) readNWISdv(y,parameterCd,
                                                            "1970-01-01","2020-12-31"))

# Rename elements with site numbers
names(rawDailyData) <- site_numbers

# Remove single site at which there is no discharge data - "03293500"
# Otherwise the following step throws an error
rawDailyData_ed <- rawDailyData[-79]

# Extract relevant information
DailyQ <- lapply(rawDailyData_ed, function(z) return(z[,c("site_no","Date","X_00060_00003","X_00060_00003_cd")]))

# Remove provisional data
DailyQ_clean <- lapply(DailyQ, function(x) return(x[which(x$X_00060_00003_cd %in% c("A","A e")),]))

################################################
## Calculate 2 year flood recurrence interval
##############################################

# testing at a single site since I'm getting some errors applying the function to all 207 sites
#Grouping maximum average daily discharge for each year
tmax.year <- DailyQ_clean$`0112400` %>% 
  mutate(Year = year(Date))

tmax.by.year <- tmax.year %>%
  group_by(Year) %>% 
  summarize(amount=max(X_00060_00003, na.rm = TRUE)) %>%
  ungroup()

#Recording the maximum discharges by year and removing N.A. values
tmaximas<-tmax.by.year$amount
tmaximas<-tmaximas[!is.na(tmaximas)]

#Sorting maxima by decreasing order
tsorted.maximas<-as.data.frame(sort(tmaximas,decreasing=T))
tsorted.maximas$rank <- seq(length=nrow(tsorted.maximas))
colnames(tsorted.maximas) <- c("Q_max","rank")

#Fit relationship
tsorted.maximas$ln_Q_max <- log(tsorted.maximas$Q_max)
tsorted.maximas$exceed_prob <- tsorted.maximas$rank/(length(tsorted.maximas$rank)+1)

#visualize
plot <- ggplot(tsorted.maximas, aes(ln_Q_max, exceed_prob))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_y_reverse()

print(plot)

#remove -Inf values
tsorted.maximas[tsorted.maximas == "-Inf"] <- NA

#extract coefficients and calc 
tmod <- lm(exceed_prob ~ ln_Q_max, data = tsorted.maximas) # There's an infinite value in there causing problems, so adding the line of code above
tl <- as.data.frame(t(as.matrix(coefficients(tmod))))
colnames(tl) <- c("int","slope")
tyr_2 <- exp((0.5 - tl$int)/tl$slope) # Adjust depending on timeframe of flood examined.
# was previously 0.5 for 2 year flood, 0.1 for 10 year flood

# Now, to build the function.
two_year_flood <- function(data){
  
  #Grouping maximum average daily discharge for each year
  max.year <- data %>% 
    mutate(Year = year(Date))
  
  max.by.year <- max.year %>%
    group_by(Year) %>% 
    summarize(amount=max(X_00060_00003, na.rm = TRUE)) %>%
    ungroup()
  
  #Recording the maximum discharges by year and removing N.A. values
  maximas<-max.by.year$amount
  maximas<-maximas[!is.na(maximas)]
  
  #Sorting maxima by decreasing order
  sorted.maximas<-as.data.frame(sort(maximas,decreasing=T))
  sorted.maximas$rank <- seq(length=nrow(sorted.maximas))
  colnames(sorted.maximas) <- c("Q_max","rank")
  
  #Fit relationship
  sorted.maximas$ln_Q_max <- log(sorted.maximas$Q_max)
  sorted.maximas$exceed_prob <- sorted.maximas$rank/(length(sorted.maximas$rank)+1)
  
  #visualize
  p <- ggplot(sorted.maximas, aes(ln_Q_max, exceed_prob))+
    geom_point()+
    geom_smooth(method = "lm")+
    scale_y_reverse()
  
  print(p)
  
  #remove NaN and Inf values
  sorted.maximas[sorted.maximas == "-Inf"] <- NA
  
  #extract coefficients and calc 
  mod <- lm(exceed_prob ~ ln_Q_max, data = sorted.maximas)
  l <- as.data.frame(t(as.matrix(coefficients(mod))))
  colnames(l) <- c("int","slope")
  yr_2 <- exp((0.5 - l$int)/l$slope) # Adjust depending on timeframe of flood examined.
  # was previously 0.5 for 2 year flood, 0.1 for 10 year flood
  
  return(yr_2)
  
}

RI_two <- ldply(lapply(DailyQ_clean, function(y) two_year_flood(y)), data.frame)
colnames(RI_two) <- c(".id","RI_2yr_Q")
RI_two$site_name <- paste("nwis_",RI_two$.id, sep = "")
RI_two <- RI_two[,c("site_name","RI_2yr_Q")]
RI_two$RI_2yr_Q_cms <- RI_two$RI_2yr_Q/35.314666212661 # convert cfs to cms

## Export - save to data folder
write_csv(RI_two, "data_working/RI_2yr_flood_182riv.csv") ## in cms

# End of script.
