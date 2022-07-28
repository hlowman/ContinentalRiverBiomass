## 207 rivers data source
## StreamLight Data Acquisition
## July 28, 2022
## Heili Lowman

# Working through third step of Stream Light data download workflow described at:
# https://psavoy.github.io/StreamLight/articles/3%20Using%20stream_light.html

# Load packages
lapply(c("tidyverse", "here",
         "StreamLightUtils", "StreamLight"), require, character.only=T)

# Preparing a driver file.

# Load in the MODIS and NALDS pre-processed files.
NLDAS_processed <- readRDS("data_raw/MODIS_processed_072822.rds")

# Stopped here because no MODIS data available at these sites.

# End of script.