## Recovery of Stream Productivity following Disturbance
## Originally created: May 11, 2023
## Heili Lowman

#### READ ME ####

# The following set of scripts will walk through the steps necessary to
# prep and send data to Beartooth as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_XX" folders have been ignored using git.ignore.
# Links to the raw data sets are provided in the 01_NWIS_RiverSelection.R file.

# If you are accessing the code via GitHub, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running the 
# code below.

# The following script will combine all data for use in post-hoc modeling
# and figure creation.

#### Setup ####

# Load necessary packages.
lapply(c("tidyverse", "lubridate", "data.table",
         "here"), require, character.only=T)

# Load necessary datasets.

# Load site-level info (Blaszczak and Appling datasets).
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
site <- fread("data_raw/site_data.tsv")

# Load in original list of data fed into the model.
dat_in <- readRDS("data_working/list_181sites_Qmaxnorm_SavoySL.rds")

# Load in model diagnostics for site-level parameters.
dat_diag <- readRDS("data_working/beartooth_181rivers_model_diags_050923.rds")

# Load in list containing all iterations of site-level parameters.
dat_out <- readRDS("data_working/beartooth_181rivers_model_params_all_iterations_050823.rds")

# Load filtered maximum accrual (yield) dataset.
dat_amax <- readRDS("data_working/maxalgalyield_152sites_051023.rds")

# Load filtered Qc:q2yr dataset.
dat_Qc <- readRDS("data_working/QcQ2_filtered_138sites_051023.rds")

# Load in nutrient data downloaded from USGS.
dat_nuts <- readRDS("data_working/USGS_WQP_nuts_aggsite_050923.rds")

# And the dataset with all HUC delineations.
site_HUC <- readRDS("data_working/NWIS_Info_181riv_HUC2_df_050923.rds")

# And the dataset with Qc exceedances/year.
dat_exc <- readRDS("data_working/Qc_exceedances_181sites_051123.rds")

#### Calculate site statistics ####

# Take list containing all input data and make into a df.
dat_in_df <- map_df(dat_in, ~as.data.frame(.x), .id="site_name")

# Calculate coefficient of variation in discharge at every site,
# as well as mean daily light availability and mean daily GPP.
dat_in_daily_means <- dat_in_df %>%
  group_by(site_name) %>%
  summarize(cvQ = (sd(Q, na.rm = TRUE)/mean(Q, na.rm = TRUE)),
            meanLight = mean(PAR_surface, na.rm = TRUE),
            meanTemp = mean(temp, na.rm = TRUE),
            meanGPP = mean(GPP, na.rm = TRUE)) %>%
  ungroup()

# Also would like to add site characteristics to this dataset 
# for plotting and linear modeling purposes.
# from Blaszczak dataset:
dat_site_info <- site_info %>%
  mutate(Order = factor(NHD_STREAMORDE)) %>%
  dplyr::select(SiteID, Lat_WGS84, Lon_WGS84, Order, NHD_AREASQKM, LU_category,
                NHD_RdDensCat, NHD_RdDensWs, NHD_PctImp2011Cat, NHD_PctImp2011Ws)
# from Appling dataset:
dat_site <- site %>%
  mutate(Canal = factor(struct.canal_flag),
         Dam = factor(struct.dam_flag)) %>%
  dplyr::select(site_name, dvqcoefs.a, dvqcoefs.b, Canal, Dam)

# Also, need to calculate stream width.
# Pull out coefficients.
coeff_ab <- dat_site %>%
  dplyr::select(site_name, dvqcoefs.a, dvqcoefs.b)
# Pull out daily discharge data for input dataset.
q_daily <- dat_in_df %>%
  dplyr::select(site_name, date, Q)
# Calculate stream width every day based on: width = a*(Q^b)
q_coeffs <- left_join(q_daily, coeff_ab)
q_coeffs$width <- q_coeffs$dvqcoefs.a*(q_coeffs$Q^q_coeffs$dvqcoefs.b)
# And summarize by site (use median).
med_width <- q_coeffs %>%
  group_by(site_name) %>%
  summarize(width_med = median(width)) %>% # in meters
  ungroup()

# And finally, edit nutrient data.
# Fix naming schema.
dat_nuts$site_name <- str_replace_all(dat_nuts$MonitoringLocationIdentifier, 'USGS-', 'nwis_')

# Widen dataset.
dat_nuts_w <- dat_nuts %>%
  pivot_wider(id_cols = "site_name",
              names_from = "CharacteristicName",
              values_from = "mean_mg_L")

# Edit and trim HUC dataset.
site_HUC_trim <- site_HUC %>%
  dplyr::select(agency_cd, site_no, huc_2) %>%
  mutate(site_name = paste(agency_cd, site_no, sep="-"))
site_HUC_trim$site_name <- str_replace_all(site_HUC_trim$site_name, 'USGS-', 'nwis_')

# Join all ancillary data together.

dat_join <- dat_in_daily_means # no longer using summer mean values
dat_join2 <- left_join(dat_join, dat_site_info, by = c("site_name" = "SiteID"))
dat_join3 <- left_join(dat_join2, dat_site)
dat_join4 <- left_join(dat_join3, med_width)
dat_join5 <- left_join(dat_join4, dat_nuts_w)
dat_join6 <- left_join(dat_join5, site_HUC_trim)
dat_join7 <- left_join(dat_join6, dat_exc)

# Export covariate data for all sites.
saveRDS(dat_join7, "data_working/covariate_data_181sites_070523.rds")

#### Join datasets for amax and Qc:Q2 separately ####

# Yield/maximum accrual (n = 152 sites)
dat_join_amax <- left_join(dat_amax, dat_join7)

# Export dataset.
saveRDS(dat_join_amax, "data_working/amax_covariates_152sites_070523.rds")

# Discharge threshold (Qc) (n = 138 sites)
dat_join_Qc <- left_join(dat_Qc, dat_join7)

# Export dataset.
saveRDS(dat_join_Qc, "data_working/Qc_covariates_138sites_070523.rds")

# End of script.
