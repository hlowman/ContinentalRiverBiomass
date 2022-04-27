## Fitting models to data
## Step FOUR in Metabolism Modeling Workflow
## April 27, 2022
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to fit the Ricker model to 43 sites for my JASM poster presentation.

# This script is actually the one run on Teton. And it is being run in parallel
# with the larger dataset to ensure I have results to put on my poster.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

#### Selecting sites of interest ####
# check to see which sites performed best in the previous model format
# divs207 <- readRDS("data_working/teton_207rivers_model_divergences_bothmodels_021322.rds")
# divs207wp <- divs207 %>% filter(model == "with P")
# divs0 <- divs207wp %>% filter(divergences == 0)
# zero_div_sites <- as.data.frame(divs0$site_name) %>% rename(site_name = `divs0$site_name`)

# load in original dataset to assign index numbers to locate within the matrices
# df <- readRDS("data_working/df_206sites_10yrQnorm.rds")
# site_name <- names(df)
# index <- seq(1, 206, by = 1)
# df2 <- as.data.frame(cbind(index, site_name))

# and now pull out indices of interest matching well-performing sites
# df_match <- left_join(zero_div_sites, df2)
# sites43 <- df_match$index
# Meaning: 7:9, 12, 15, 16, 22, 23, 25, 32, 33, 39, 49, 55, 60, 66, 70:72, 76, 79, 81, 98, 101, 115, 116, 119, 120, 132:134, 140, 142, 157, 163:166, 169, 172, 176, 189, 200

# Load in matrices with all 206 sites
# light_mx <- readRDS("data_working/light_mx206.rds")
# gpp_mx <- readRDS("data_working/gpp_mx206.rds")
# gpp_sd_mx <- readRDS("data_working/gpp_sd_mx206.rds")
# tQ_mx <- readRDS("data_working/tQ_mx206.rds")
# e_mx <- readRDS("data_working/e_mx206.rds")

# Subset to 43 sites using indices above.
# light_mx43 <- light_mx[1:2926, c(7:9, 12, 15, 16, 22, 23, 25, 32, 33, 39, 49, 55, 60, 66, 70:72, 76, 79, 81, 98, 101, 115, 116, 119, 120, 132:134, 140, 142, 157, 163:166, 169, 172, 176, 189, 200)]
# gpp_mx43 <- gpp_mx[1:2926, c(7:9, 12, 15, 16, 22, 23, 25, 32, 33, 39, 49, 55, 60, 66, 70:72, 76, 79, 81, 98, 101, 115, 116, 119, 120, 132:134, 140, 142, 157, 163:166, 169, 172, 176, 189, 200)]
# gpp_sd_mx43 <- gpp_sd_mx[1:2926, c(7:9, 12, 15, 16, 22, 23, 25, 32, 33, 39, 49, 55, 60, 66, 70:72, 76, 79, 81, 98, 101, 115, 116, 119, 120, 132:134, 140, 142, 157, 163:166, 169, 172, 176, 189, 200)]
# tQ_mx43 <- tQ_mx[1:2926, c(7:9, 12, 15, 16, 22, 23, 25, 32, 33, 39, 49, 55, 60, 66, 70:72, 76, 79, 81, 98, 101, 115, 116, 119, 120, 132:134, 140, 142, 157, 163:166, 169, 172, 176, 189, 200)]
# e_mx43 <- e_mx[1:2926, c(7:9, 12, 15, 16, 22, 23, 25, 32, 33, 39, 49, 55, 60, 66, 70:72, 76, 79, 81, 98, 101, 115, 116, 119, 120, 132:134, 140, 142, 157, 163:166, 169, 172, 176, 189, 200)]

# Cross-checked a site on each data type to be sure the correct data was being pulled.

# Export these matrices for use in the actual Teton run.
# saveRDS(light_mx43, "data_working/light_mx43.rds")
# saveRDS(gpp_mx43, "data_working/gpp_mx43.rds")
# saveRDS(gpp_sd_mx43, "data_working/gpp_sd_mx43.rds")
# saveRDS(tQ_mx43, "data_working/tQ_mx43.rds")
# saveRDS(e_mx43, "data_working/e_mx43.rds")

## Source data - all in matrix format for new random effects structure
line_lengths <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/line_lengths206.rds")
light_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/light_mx43.rds")
gpp_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/gpp_mx43.rds")
gpp_sd_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/gpp_sd_mx43.rds")
tQ_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/tQ_mx43.rds")
e_mx <- readRDS("/project/modelscape/users/hlowman/jobscripts/teton_jasm/e_mx43.rds")


####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data for partial run of 43 sites (max of 2926 observations)
stan_data_l <- list(sites = 43, # number of sites 
                    Nobs = 2926, # max number of observations (days)
                    Ndays = line_lengths[c(7:9, 12, 15, 16, 22, 23, 25, 32, 33, 39, 49, 55, 60, 66, 70:72, 76, 79, 81, 98, 101, 115, 116, 119, 120, 132:134, 140, 142, 157, 163:166, 169, 172, 176, 189, 200),], # number of observations per site
                    light = light_mx, # standardized light data
                    GPP = gpp_mx, # standardized GPP estimates
                    GPP_sd = gpp_sd_mx, # standardized GPP standard deviations
                    tQ = tQ_mx, # 10 yr flood standardized discharge
                    new_e = e_mx) # indices denoting when to reinitialize biomass estimation

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model

# sets initial values of c and s to help chains converge
init_Ricker <- function(...) {
  list(csite = rep(0.5,length.out = 43), 
       ssite = rep(0.5,length.out = 43),
       rsite = rep(0.5,length.out = 43),
       lsite = rep(0.5,length.out = 43)) # new values as of apr 2022
}

## export results

PM_outputlist_Ricker <- stan("/project/modelscape/users/hlowman/jobscripts/teton_jasm/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_wP_re43.stan",
                             data = stan_data_l, 
                             chains = 3,
                             iter = 5000,
                             init = init_Ricker,
                             control = list(max_treedepth = 12))

saveRDS(PM_outputlist_Ricker, "/project/modelscape/users/hlowman/jobresults/teton_jasm/teton_43rivers_output_Ricker_2022_04_27.rds")

# End of script.
