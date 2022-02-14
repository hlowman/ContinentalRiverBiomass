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

# saveRDS(PM_outputlist_Ricker_wP, "data_working/supp_figures/potomac_output_Ricker_wP_2022_01_31.rds")
PM_outputlist_Ricker_wP <- readRDS("data_working/supp_figures/potomac_output_Ricker_wP_2022_01_31.rds")

# without P term
PM_outputlist_Ricker_noP <- lapply(stan_data_l,
                                   function(x) stan("code/supp_figures/Stan_ProductivityModel2_Ricker_fixedinit_obserr_ts_noP_02_2022.stan",
                                                data = x, chains = 3,iter = 5000,
                                                init = init_Ricker,
                                                control = list(max_treedepth = 12)))

# saveRDS(PM_outputlist_Ricker_noP, "data_working/supp_figures/potomac_output_Ricker_noP_2022_01_31.rds")
PM_outputlist_Ricker_noP <- readRDS("data_working/supp_figures/potomac_output_Ricker_noP_2022_01_31.rds")

## modify results
# Going to create functions to pull out data of interest.
extract_params <- function(df){
  extract(df, c("r", "lambda", "s", "c", "sig_p", "sig_o"))
}

extract_params_noP <- function(df){
  extract(df, c("r", "lambda", "sig_p", "sig_o"))
}

# With P
# And now map this to the entire output list.
data_out_params <- map(PM_outputlist_Ricker_wP, extract_params)

# And create a dataframe
params_wP_df <- map_df(data_out_params, ~as.data.frame(.x), .id="site_name") %>%
  mutate(model = "with P") %>%
  mutate(k = (-1*r)/lambda)

# And summarize
params_means_wP <- params_wP_df %>%
  group_by(site_name) %>%
  summarize(r_mean = mean(r), 
            l_mean = mean(lambda), 
            s_mean = mean(s), 
            c_mean = mean(c), 
            sigp_mean = mean(sig_p), 
            sigo_mean = mean(sig_o))

# Without P
data_out_params_noP <- map(PM_outputlist_Ricker_noP, extract_params_noP)

params_noP_df <- map_df(data_out_params_noP, ~as.data.frame(.x), .id="site_name") %>%
  mutate(model = "without P") %>%
  mutate(k = (-1*r)/lambda)

params_means_noP <- params_noP_df %>%
  group_by(site_name) %>%
  summarize(r_mean = mean(r), 
            l_mean = mean(lambda), 
            sigp_mean = mean(sig_p), 
            sigo_mean = mean(sig_o))

#### GPP Predictions ####

# Using code from Joanna's script "Predicted_ProductivityModel_Ricker.R".

## Growth Model 3 - Data simulation

# Function to predict GPP with P term
PM_Ricker_wP <- function(r, lambda, s, c, sig_p, sig_o, df) {
  
  ## Data
  Ndays <- length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value
  
  ## Vectors for model output of P, B, pred_GPP
  P <- numeric(Ndays)
  P[1] <- 1
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*100*(tQ[i] - c)))
  }
  
  B <- numeric(Ndays)
  B[1] <- log(GPP[1]/light[1])
  
  pred_GPP <- numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] <- MCMCglmm::rtnorm(1, 
                             mean = (B[j-1] + r + lambda*exp(B[j-1]))*P[j], 
                             sd = sig_p, upper = 5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, 
                                    mean = light[i]*exp(B[i]), 
                                    sd = sig_o, lower=0.01)
  }
  
  return(pred_GPP)
}

# Function to predict GPP withOUT P term
PM_Ricker_woP <- function(r, lambda, sig_p, sig_o, df) {
  
  ## Data
  Ndays <- length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  # tQ <- df$tQ # discharge standardized to max value
  
  ## Vectors for model output of B, pred_GPP
  B <- numeric(Ndays)
  B[1] <- log(GPP[1]/light[1])
  
  pred_GPP <- numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] <- MCMCglmm::rtnorm(1, 
                             mean = (B[j-1] + r + lambda*exp(B[j-1])), 
                             sd = sig_p, upper = 5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, 
                                    mean = light[i]*exp(B[i]), 
                                    sd = sig_o, lower=0.01)
  }
  
  return(pred_GPP)
}

# Predict GPP for 4 cases:
# (1) Potomac with P
# values from params_means_wP
predPOT_P <- PM_Ricker_wP(r = 0.307082, lambda = -0.02712014, s = 1.74471630, 
                          c = 0.1489799, sig_p = 0.2092323, sig_o = 0.90304684, 
                          df = l$nwis_01608500)

# (2) Potomac without P
# values from params_means_noP
predPOT_noP <- PM_Ricker_woP(r = 0.1692938, lambda = -0.01653095, 
                             sig_p = 0.2072508, sig_o = 1.01722985, 
                             df = l$nwis_01608500)

# (3) Reedy with P
# values from params_means_wP
predREED_P <- PM_Ricker_wP(r = -1.662759, lambda = -0.98981501, s = 0.03080441, 
                           c = 0.3646432, sig_p = 2.0452212, sig_o = 0.04770292, 
                           df = l$nwis_02266300)

# (4) Reedy without P
# values from params_means_noP
predREED_noP <- PM_Ricker_woP(r = 0.3703312, lambda = -1.24159157,
                              sig_p = 2.1290656, sig_o = 0.05169693, 
                              df = l$nwis_02266300)

# NOTE TO FUTURE SELF - need to add confidence intervals to this figure
# using Joanna's code at "Biomass2_WSpredictions.R" or "Biomass2a_Fig.R"

rivers_gpp <- rivers %>%
  select(date, site_name, GPP)

rivers_dates <- rivers_gpp$date

# Bind predictions from above with dates and sites
# with P term
predPOT_P_df <- as.data.frame(predPOT_P) %>%
  mutate(site_name = "nwis_01608500") %>%
  rename(predGPP_P = predPOT_P)

predREED_P_df <- as.data.frame(predREED_P) %>%
  mutate(site_name = "nwis_02266300") %>%
  rename(predGPP_P = predREED_P)

pred_P_df <- rbind(predPOT_P_df, predREED_P_df)
pred_P_df <- cbind(rivers_dates, pred_P_df)

# without P term
predPOT_noP_df <- as.data.frame(predPOT_noP) %>%
  mutate(site_name = "nwis_01608500") %>%
  rename(predGPP_noP = predPOT_noP)

predREED_noP_df <- as.data.frame(predREED_noP) %>%
  mutate(site_name = "nwis_02266300") %>%
  rename(predGPP_noP = predREED_noP)

pred_noP_df <- rbind(predPOT_noP_df, predREED_noP_df)
pred_noP_df <- cbind(rivers_dates, pred_noP_df)

gpp_with_predictions1 <- full_join(pred_P_df, pred_noP_df) %>%
  rename(date = rivers_dates)
gpp_with_predictions <- full_join(gpp_with_predictions1, rivers_gpp) %>%
  mutate(Date = ymd(date)) %>%
  select(site_name, Date, GPP, predGPP_P, predGPP_noP) %>%
  pivot_longer(cols = starts_with("pred"), names_to = "model", values_to = "predGPP") %>%
  mutate(model_f = factor(model, levels = c("predGPP_P", "predGPP_noP")))

# new facet labels for sites
site.labs <- c("South Branch Potomac River (WV)", "Reedy Creek (FL)")
names(site.labs) <- c("nwis_01608500", "nwis_02266300")

model.labs <- c("With P Term", "Without P Term")
names(model.labs) <- c("predGPP_P", "predGPP_noP")

supp_fig <- ggplot(gpp_with_predictions)+
  geom_point(aes(x = Date, y = GPP), color="black", size=1)+
  geom_line(aes(x = Date, y = predGPP, color=model_f), alpha=0.8, size=1)+
  scale_color_manual(values = c("#3793EC", "#6CA184"))+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x = "Date") +
  scale_x_date(date_labels = "%b") +
  facet_grid(site_name~model_f, scales = "free",
             labeller = labeller(site_name = site.labs, model_f = model.labs)) +
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
#        height = 12,
#        units = "cm"
# )

# End of script.
