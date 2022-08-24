# code for explanatory variable exploration

library("tidyverse")
library("here")
library("GGally")
library("glmmTMB")
library("MuMIn")
library("effects")
library("DHARMa")

# load in full hypoxia site information dataset
hypoxia_dataset <- read_csv("data_raw/GRDO_GEE_HA_NHD.csv")

# trim down to sites of interest and select only columns of interest
hypoxia_PC_trim <- hypoxia_dataset %>%
  filter(DB_Source == "PC") %>%
  select(SiteID, Lat_WGS84, Lon_WGS84, ele_mt_cav, CATCH_SKM, NHD_AREASQKM,
         dis_m3_pyr, pre_mm_cyr, glc_cl_cmj, LU_category, NHD_RdDensCat, tmp_dc_cyr)

#### Correlation ####

# check correlation of variables independent of rmax
(corr_sites <- ggpairs(hypoxia_PC_trim %>% select(-c(SiteID, Lat_WGS84, Lon_WGS84,
                                                  ele_mt_cav, NHD_AREASQKM, glc_cl_cmj))))

# air temperature and precip seem correlated, as do watershed size and precip (although I'm not sure as to why for the second one)

#### Linear Model ####

# log-transform prior to running model

m1 <- glmmTMB(rmax ~ CATCH_SKM + LU_category + # subsidy variables
                dis_m3_pyr, # + mean annual PAR #disturbance variables
              data = my_data_here)

m2 <- glmmTMB(rmax ~ LU_category + # subsidy variables
                dis_m3_pyr + pre_mm_cyr, # + mean annual PAR #disturbance variables
              data = my_data_here)

m3 <- glmmTMB(rmax ~ CATCH_SKM + LU_category + # subsidy variables
                dis_m3_pyr + tmp_dc_cyr, # + mean annual PAR #disturbance variables
              data = my_data_here)

AICc(m1, m2, m3)

# Output of the model.
# Note, summary() function looks at contrasts between singular effects.
summary(m1)

# Checking residuals.
plot(m1, col=1)
qqnorm(m1)

# Final results.
# anova() function looks at contrasts across all effects.
anova(m1)

# Post-hoc:
mHSD <- glht(m1, linfct=mcp(LU_category="Tukey")) # Run a Tukey's post hoc analysis on land use.
summary(mHSD) 

# End of script.
