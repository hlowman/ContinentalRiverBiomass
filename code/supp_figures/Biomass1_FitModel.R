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

reedy <- df$nwis_02266300 %>%
  filter(year == 2012)

rivers <- rbind(potomac, reedy)

## split list by ID
l <- split(rivers, rivers$site_name)

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

stan_data_l <- lapply(l, function(x) stan_data_compile(x))

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
PM_outputlist_Ricker_wP <- lapply(stan_data_l,
                                  function(x) stan("code/supp_figures/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_02_2022.stan",
                                                data = x, chains = 3,iter = 5000,
                                                init = init_Ricker,
                                                control = list(max_treedepth = 12)))

saveRDS(PM_outputlist_Ricker_wP, "data_working/supp_figures/potomac_output_Ricker_wP_2022_01_31.rds")

# without P term
PM_outputlist_Ricker_noP <- lapply(stan_data_l,
                                   function(x) stan("code/supp_figures/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_noP_02_2022.stan",
                                                data = x, chains = 3,iter = 5000,
                                                init = init_Ricker,
                                                control = list(max_treedepth = 12)))

saveRDS(PM_outputlist_Ricker_noP, "data_working/supp_figures/potomac_output_Ricker_noP_2022_01_31.rds")

## modify results
# Going to create a function of the above to pull out data of interest from
# all sites.
extract_params <- function(df){
  extract(df, c("pred_GPP"))
}

# With P
# And now map this to the entire output list.
data_out_params <- map(PM_outputlist_Ricker_wP, extract_params)
# the above line of code sometimes doesn't play nicely if R has been up and running
# for awhile, so the fix is to exit RStudio and reopen the project/file

# And create a dataframe
params_wP_df <- map_df(data_out_params, ~as.data.frame(.x), .id="site_name") %>%
  mutate(model = "with P")

# Without P
data_out_params_noP <- map(PM_outputlist_Ricker_noP, extract_params)

params_noP_df <- map_df(data_out_params_noP, ~as.data.frame(.x), .id="site_name") %>%
  mutate(model = "without P")

# NOTE TO FUTURE SELF - need to add intervals to this figure

rivers_gpp <- rivers %>%
  select(date, site_name, GPP)

rivers_dates <- rivers_gpp$date

# Bind model outputs from above together
params_all <- rbind(params_wP_df, params_noP_df)

params_gpp <- params_all %>%
  group_by(model, site_name) %>% # group by site and model run
  summarise_all(list(mean)) %>% # calculate means
  #select(c(model, contains('GPP'))) %>% # select model column and all GPP columns
  pivot_longer(cols = pred_GPP.1:pred_GPP.352, names_to = "day", values_to = "pred_GPP") %>% # pivot longer
  drop_na(pred_GPP) %>% # remove days on which there were no predictions
  ungroup()

params_gpp <- params_gpp %>%
  mutate(date = rep(rivers_dates, 2))  # add dates back in

# new facet labels for sites
site.labs <- c("South Branch Potomac River (WV)", "Reedy Creek (FL)")
names(site.labs) <- c("nwis_01608500", "nwis_02266300")

supp_fig <- ggplot()+
  geom_point(data = rivers_gpp, aes(x = date, y = GPP), color="gray60", size=2)+
  geom_line(data = params_gpp , aes(x = date, y = pred_GPP, color = model), size=1)+
  scale_color_manual(values = c("#3793EC", "#6CA184"))+
  #geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="darkolivegreen4")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x = "Date")+
  facet_grid(model~site_name, scales = "free",
             labeller = labeller(site_name = site.labs)) +
  theme(strip.background = element_rect(fill = NA),
        panel.background = element_rect(color = "black", fill=NA, size=1),
        legend.position = "none",
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=12))

supp_fig

# ggsave(("figures/supp_figures/nwis_01608500_02266300_with_without_P.png"),
#        width = 20,
#        height = 20,
#        units = "cm"
# )

# End of script.
