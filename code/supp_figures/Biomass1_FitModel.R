## Fitting models to data
## Step FOUR in Metabolism Modeling Workflow
## January 31, 2022
## Heili Lowman

# I'll be modifying Joanna's code from the RiverBiomass repository
# to fit the Ricker model to 1 site of data to generate results for
# a supplementary figure demonstrating the utility of removing P.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data
df <- readRDS("data_working/df_207sites.rds")

potomac <- df$nwis_01608500 %>%
  filter(year == 2012)

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data
stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), 
               light = x$light_rel, 
               GPP = x$GPP,
               GPP_sd = x$GPP_sd, 
               tQ = x$tQ,
               new_e = x$new_e) # new column for reinitialization of B[j] values
               
  return(data)
}

stan_data_l <- stan_data_compile(potomac)

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

# Latent Biomass (Ricker population) Model

# sets initial values of c and s to help chains converge
init_Ricker <- function(...) {
  list(c = 0.5, s = 0.5) # new values as of jan 2022
}

## export results
# with P term
PM_outputlist_Ricker_wP <- stan("code/supp_figures/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_02_2022.stan",
                                                data = stan_data_l, chains = 3,iter = 5000,
                                                init = init_Ricker,
                                                control = list(max_treedepth = 12))

saveRDS(PM_outputlist_Ricker_wP, "data_working/supp_figures/potomac_output_Ricker_wP_2022_01_31.rds")

# without P term
PM_outputlist_Ricker_noP <- stan("code/supp_figures/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_noP_02_2022.stan",
                                data = stan_data_l, chains = 3,iter = 5000,
                                init = init_Ricker,
                                control = list(max_treedepth = 12))

saveRDS(PM_outputlist_Ricker_noP, "data_working/supp_figures/potomac_output_Ricker_noP_2022_01_31.rds")

#### Output processing and figure creation ####
params_wP <- extract(PM_outputlist_Ricker_wP, c("r","lambda",
                                                 "B","pred_GPP","sig_p","sig_o"))
# not pulling out s, c, and P so dataframes match in orientation for joining below

# And create a dataframe
params_wP_df <- as.data.frame(params_wP) %>%
  # and add "K" to it, calculating for each individual iteration
  mutate(k = (-1*r)/lambda) %>%
  mutate(model = "with P")

params_noP <- extract(PM_outputlist_Ricker_noP, c("r","lambda",
                                                "B","pred_GPP","sig_p","sig_o"))

# And create a dataframe
params_noP_df <- as.data.frame(params_noP) %>%
  # and add "K" to it, calculating for each individual iteration
  mutate(k = (-1*r)/lambda) %>%
  mutate(model = "without P")

# NOTE TO FUTURE SELF - need to add intervals to this figure

potomac_gpp <- potomac %>%
  select(date, GPP)

dates <- potomac_gpp$date

# Bind model outputs from above together
params_all <- rbind(params_wP_df, params_noP_df)

params_gpp <- params_all %>%
  group_by(model) %>% # group by model run
  summarise_all(list(mean)) %>% # calculate means
  select(c(model, contains('GPP'))) %>% # select model column and all GPP columns
  pivot_longer(!model, names_to = "day", values_to = "pred_GPP") %>% # pivot longer
  mutate(date = rep(dates, 2)) # add dates back in

supp_fig <- ggplot()+
  geom_line(data = params_gpp, aes(x = date, y = pred_GPP, color = model), size=1)+
  scale_color_manual(values = c("#4CA49E", "#6B6D9F"))+
  geom_point(data = potomac_gpp, aes(x = date, y = GPP), color="black", size=2)+
  #geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="darkolivegreen4")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x = "Date")+
  theme(panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12))

supp_fig

# ggsave(("figures/supp_figures/nwis_01608500_with_without_P.png"),
#        width = 20,
#        height = 10,
#        units = "cm"
# )

# End of script.
