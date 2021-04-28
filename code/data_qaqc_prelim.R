# Data QAQC
# April 28, 2021
# Heili Lowman

# This script will be use for preliminary filtering
# of the Appling et al. 2018 dataset (citation below)
# to create a dataset of usable metabolism estimates
# based on the following criteria:

# (1) Sites should have a minimum of 275 days of data/year.

# (2) Sites should have data gaps no longer than 14 days.

# (3) Rhat values for GPP estimates should be less than 1.05.


# Setup -------------------------------------------------------------------

# Load packages.
library(tidyverse)
library(lubridate)
library(here)

# Double check working directory
here()

# Load datasets.
dailydat <- read_tsv(here("data_raw", 
                          "daily_predictions.tsv"))

