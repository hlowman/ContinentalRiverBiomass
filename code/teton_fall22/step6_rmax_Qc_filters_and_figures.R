## Resilience of Stream Productivity to Disturbance
## October 13, 2022
## Heili Lowman

# The following set of scripts will walk through the steps necessary to
# prep and send data to Teton as well as process the model outputs.

# Much of this code has been modified from the RiverBiomass repository
# found at: https://github.com/jrblaszczak/RiverBiomass 

# Please note, the "data_raw" and "data_working" folders have been ignored
# using git.ignore, so links to the raw data sets are provided in the step1
# file. If you are accessing the code via GitHub, these will need to be 
# downloaded and added to a folder of the appropriate name prior to running
# the code.

#### Setup ####

# Load necessary packages.
lapply(c("tidyverse", "lubridate", "data.table",
         "rstan","bayesplot","shinystan", "here",
         "GGally", "glmmTMB", "MuMIn", "effects",
         "DHARMa", "lme4", "multcomp", "patchwork",
         "calecopal", "viridis"), require, character.only=T)

# Load necessary datasets.
# Load site-level info (hypoxia and Appling datasets).
site_info <- read_csv("data_raw/GRDO_GEE_HA_NHD_2021_02_07.csv")
site <- fread("data_raw/site_data.tsv")

# Load in original list of data fed into the model.
dat_in <- readRDS("data_working/list_182sites_Qmaxnorm_allSL.rds")

# Load in model diagnostics for site-level parameters.
dat_diag <- readRDS("data_working/teton_182rivers_model_diags_101522.rds")

# Load in list containing all iterations of site-level parameters.
dat_out <- readRDS("data_working/teton_182rivers_model_params_all_iterations_101522.rds")

# Load in 2 year flood values.
dat_2yr <- read_csv("data_working/RI_2yr_flood_182riv.csv")

# Load in nutrient data downloaded from USGS.
dat_nuts <- readRDS("data_working/USGS_WQP_nuts_aggsite_110222.rds")

#### Data Prep ####

# Take list containing all input data and make into a df.
dat_in_df <- map_df(dat_in, ~as.data.frame(.x), .id="site_name")

# Take list containing all iterations of parameters and make into a df.
dat_out_df <- map_df(dat_out, ~as.data.frame(.x), .id="site_name")

#### Value filter for rmax ####

# Negative rmax values are not biologically reasonable, so I've 
# removed them.

# First, calculate mean rmax values at all the sites.
# dat_out_rmean <- dat_out_df %>%
#   group_by(site_name) %>%
#   summarize(r_mean = mean(r)) %>%
#   ungroup()

# Instead of calculating the mean, I'll be using median rmax values.
dat_out_rmed <- dat_out_df %>%
  group_by(site_name) %>%
  summarize(r_med = median(r)) %>%
  ungroup()

# And remove negative values.
dat_out_rmed_pos <- dat_out_rmed %>%
  filter(r_med > 0) # Removes 12 sites.

#### Rhat filter for rmax ####

# Before proceeding with the first step on my analyses, I will be filtering out 
# sites at which the model did not converge well for the rmax parameter.
# Sites with Rhat > 1.05 will not pass muster.

dat_diag_rfilter1 <- dat_diag %>%
  filter(parameter == "r") %>%
  filter(Rhat < 1.05) # 12 sites drop off

#### normRMSE filter ####

# NOTE - NOT USING nRMSE FILTER FOR NOV MEETING - TOO COMPUTATIONALLY INTENSIVE

# Next, I will be examining sites that do not do a good job of predicting
# GPP using the original data.

# Using code from Joanna's scripts "Predicted_ProductivityModel_Ricker.R"
# and "Biomass2_WSpredictions.R".

# First, need to create the function for predicting GPP.
PM_Ricker <- function(r, lambda, s, c, sig_p, sig_o, df) {
  
  ## Data
  Ndays <- length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$Q_rel # discharge standardized to max value
  new_e <- df$new_e
  
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
  for(j in 2:Ndays){
    # adding in section for my re-initialization functionality
    if (new_e[j]==1) {
      
      B[j] ~ MCMCglmm::rtnorm(1, mean = log(GPP[j]/light[j])) }
    
    else {
      
      B[j] <- MCMCglmm::rtnorm(1, mean = (B[j-1] + r + lambda*exp(B[j-1]))*P[j],
                               sd = sig_p, upper = 5) }
    
  }
  
  for(i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), 
                                    sd = sig_o, lower = 0.01)
  }
  
  return(pred_GPP)
}

# Next, need to write the function with which to perform the simulation.
Ricker_sim_fxn <- function(y, x){
  # identify data
  output <- y # Teton/stan output
  df <- x # original data input
  
  # extracted parameters from STAN output already
  pars <- output
  
  # create empty matrix with days of GPP x length of iterations to receive values
  simmat <- matrix(NA, length(df$GPP), length(unlist(pars$sig_p)))
  rmsemat <- matrix(NA, length(df$GPP), 1)
  
  # simulate pred_GPP holding a parameter set for a given iteration constant
  # and then predict forward for a site's timeseries (i.e., length(df$GPP))
  for(i in 1:length(pars$r)){
    simmat[,i] <- PM_Ricker(pars$r[i], pars$lambda[i], pars$s[i], pars$c[i], pars$sig_p[i], pars$sig_p[i], df)
    rmsemat[i] <- sqrt(sum((simmat[,i] - df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat, rmsemat)
  return(l)
  
}

# And finally, apply the function to my data.
# Please note, this takes HOURS to run, so start this early in the day, and
# come back to it.
Ricker_sim_182sites <- mapply(Ricker_sim_fxn, dat_out, dat_in)

# For some reason, it's yielding two values, so let's see what's happening here.
Ricker_sim_1site <- Ricker_sim_fxn(dat_out$nwis_01124000, dat_in$nwis_01124000)
# Ok, so it's yielding pred_GPP in the first part and rmse in the second.

# Now, making a longer list to see how it spits out multiple sites to decipher
# my larger output structure.
dat_in2 <- dat_in[1:2]
dat_out2 <- dat_out[1:2]

Ricker_sim_2site <- mapply(Ricker_sim_fxn,dat_out2, dat_in2)
# Viewing this yields nothing, because it's a matrix >_<

# So, for reference:
predGPP <- Ricker_sim_1site[[1]]
rmse <- Ricker_sim_1site[[2]]

# But when this is made larger, the portions of the matrix can be accessed by
# indexing by odd and event indices.
# ODD = predGPP
# EVEN = rmse
# So, making a list of even numbers to pull out rmse values.
my_values <- seq(from = 2, to = 364, by = 2)
rmse_182sites <- Ricker_sim_182sites[my_values]

# Adding the nRMSE calculation into the function above didn't play nicely with
# the list that existed, so calculating outside instead.
nRMSE_fxn <- function(df, df_orig){
  
  # Calculate the mean RMSE value for each site.
  nRMSE <- mean(df)/(max(df_orig$GPP) - min(df_orig$GPP))
  
}

nRMSE_182sites <- mapply(nRMSE_fxn, rmse_182sites, dat_in)

nRMSE_182sitesdf <- as.data.frame(nRMSE_182sites) %>%
  mutate("site_name" = names(dat_in)) %>%
  rename("nRMSE" = "nRMSE_182sites")

# Export both sets of results.
# trying to export this first file caused the server to freeze, so only
# exported the nRMSE file for now.
#saveRDS(Ricker_sim_182sites, "data_working/Sim_Ricker_182sites_101922.rds")
saveRDS(nRMSE_182sitesdf, "data_working/nRMSE_182sites_101922.rds")

#### Additional data calculations ####

# Next, append the positive rmax values to the Rhat filter to remove
# appropriate sites.
dat_out_rmed_Rhat <- inner_join(dat_diag_rfilter1, dat_out_rmed_pos) 
# 159 sites remaining

# Not removing any sites based on RMSE values for the time being.
# REVISIT THIS following co-authors meeting in Nov.

# Finally, calculate coefficient of variation in discharge at every site,
# as well as mean daily light availability, and add to the dataset for plotting
# purposes.
dat_in_cvq_L <- dat_in_df %>%
  group_by(site_name) %>%
  summarize(cvQ = (sd(Q, na.rm = TRUE)/mean(Q, na.rm = TRUE)),
            meanL = mean(PAR_surface, na.rm = TRUE),
            meanGPP = mean(GPP, na.rm = TRUE)) %>%
  ungroup()

# Need to also calculate summer light and temperature means
dat_in_summer <- dat_in_df %>%
  filter(DOY > 164) %>% # summer solstice-ish
  filter(DOY < 259) %>% # autumn equinox-ish
  group_by(site_name) %>%
  summarize(summerL = sum(PAR_surface, na.rm = TRUE), # cumulative summer light
            summerT = mean(temp, na.rm = TRUE)) %>% # mean summer temperature
  ungroup()

# And, append this to the larger dataset.
dat_out_yas <- left_join(dat_out_rmed_Rhat, dat_in_cvq_L)
dat_out_queen <- left_join(dat_out_yas, dat_in_summer)

# Also would like to add site characteristics to this dataset for plotting.
# from hypoxia dataset:
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
  dplyr::select(site_name, Q)
# Calculate stream width every day based on: width = a*(Q^b)
q_coeffs <- full_join(q_daily, coeff_ab)
q_coeffs$width <- q_coeffs$dvqcoefs.a*(q_coeffs$Q^q_coeffs$dvqcoefs.b)
# And summarize by site (use median).
med_width <- q_coeffs %>%
  group_by(site_name) %>%
  summarize(width_med = median(width)) %>%
  ungroup() %>% # oh so this has 602 sites to reflect the orig site info dataset so,
  drop_na(width_med) # dropping NAs (two sites have no data?)

# And finally, edit nutrient data.
dat_nuts$site_name <- str_replace_all(dat_nuts$MonitoringLocationIdentifier, 'USGS-', 'nwis_')

# Widen dataset.
dat_nuts_w <- dat_nuts %>%
  pivot_wider(id_cols = "site_name",
              names_from = "CharacteristicName",
              values_from = "mean_mg_L")

# And append.
dat_out_join1 <- left_join(dat_out_queen, dat_site_info,
                          by = c("site_name" = "SiteID"))
dat_out_join2 <- left_join(dat_out_join1, dat_site)
dat_out_join3 <- left_join(dat_out_join2, med_width)
dat_out_full <- left_join(dat_out_join3, dat_nuts_w)

#### rmax Figures ####

# Distribution of rmax values:
(fig1 <- ggplot(dat_out_full, aes(x = r_med)) +
  geom_histogram(bins = 60, alpha = 0.8, 
                 fill = "#9E8ABC", color = "#9E8ABC") +
  labs(x = expression(Maximum~Growth~Rate~(r[max])),
       y = "Count") +
  theme_bw())

# Mean daily GPP vs. rmax: X axis LOG SCALED
(fig1.1 <- ggplot(dat_out_full, aes(x = meanGPP, y = r_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#A698D3") +
    scale_x_log10() +
    labs(y = expression(Maximum~Growth~Rate~(r[max])),
         x = expression(Mean~Daily~GPP~(gO[2]~m^-2~d^-1))) +
    theme_bw())

# CV of Discharge vs. rmax:
(fig2 <- ggplot(dat_out_full, aes(x = cvQ, y = r_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#9092AD") +
    labs(x = expression(CV[Q]),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

(fig2.2 <- ggplot(dat_out_full, aes(x = cvQ, y = r_med)) +
    geom_point(alpha = 0.8, size = 3) +
    labs(x = expression(CV[Q]),
         y = expression(r[max])) +
    theme_bw() +
    theme(text = element_text(family = "serif", size = 40)))

# Mean Daily Light Availability vs. rmax:
(fig3 <- ggplot(dat_out_full, aes(x = meanL, y = r_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#9E8ABC") +
    labs(x = expression(Mean~Daily~PAR~(mol~m^-2~d^-1)),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

(fig3.1 <- ggplot(dat_out_full, aes(x = summerL, y = r_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#8F8D88") +
    labs(x = expression(Cumulative~Summer~PAR~(mol~m^-2~d^-1)),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

(fig3.2 <- ggplot(dat_out_full, aes(x = summerT, y = r_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#A393CA") +
    labs(x = expression(Mean~Summer~Temperature~(Celsius)),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# Stream Order vs. rmax: Removing singular site w/o order info for now.
(fig4 <- ggplot(dat_out_full %>%
                  filter(!is.na(Order)), aes(x = Order, y = r_med)) +
    geom_boxplot(alpha = 0.6, color = "#9494B4", fill = "#9494B4") +
    labs(x = expression(Stream~Order),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# Stream Width vs. rmax: note, x axis LOG SCALED
(fig4.1 <- ggplot(dat_out_full, aes(x = width_med, y = r_med)) +
    geom_point(alpha = 0.6, size = 3, color = "#808C91") +
    scale_x_log10() + 
    labs(x = expression(Stream~Width~(m)),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# Also did a quick gut check of order vs. width
plot(as.numeric(dat_out_full$Order), dat_out_full$width_med)

# First outlier for the first order streams is Quinnipiac River, CT, which
# appears to drain straight from Hanover Pond, hence why it's so bit.
# The second outlier is Alcovy River, GA, which does NOT look like a first
# order stream; it drains a few other streams but *also* has some sort of
# retention pond upstream, which may also have reset the order classification.

# Latitude vs. rmax:
(fig5 <- ggplot(dat_out_full, aes(x = Lat_WGS84, y = r_med)) +
    geom_point(alpha = 0.6, size = 3, color = "#938E86") +
    labs(x = expression(Latitude),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# Longitude vs. rmax:
# (fig6 <- ggplot(dat_out_full, aes(x = Lon_WGS84, y = r_mean)) +
#     geom_point(alpha = 0.6, size = 3, color = "#A18F7E") +
#     labs(x = expression(Longitude),
#          y = expression(Maximum~Growth~Rate~(r[max]))) +
#     theme_bw())

# Catchment size vs. rmax: note, missing Miss. R. and x axis LOG SCALED
(fig7 <- ggplot(dat_out_full, aes(x = NHD_AREASQKM, y = r_med)) +
    geom_point(alpha = 0.6, size = 3, color = "#A99CD9") +
    scale_x_log10() +
    labs(x = expression(Watershed~Area~(km^2)),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# Land use vs. rmax:
# (fig8 <- ggplot(dat_out_full, aes(x = LU_category, y = r_mean)) +
#     geom_boxplot(alpha = 0.6, color = "#A5BA92", fill = "#A5BA92") +
#     labs(x = expression(Land~Use),
#          y = expression(Maximum~Growth~Rate~(r[max]))) +
#     theme_bw())

(fig9 <- ggplot(dat_out_full, aes(x = NHD_RdDensCat, y = r_med)) +
    geom_point(alpha = 0.6, size = 3, color = "#A6987F") +
    labs(x = expression(Road~Density~by~Catchment~(km/km^2)),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

(fig10 <- ggplot(dat_out_full, aes(x = NHD_RdDensWs, y = r_med)) +
    geom_point(alpha = 0.6, size = 3, color = "#A6987F") +
    labs(x = expression(Road~Density~by~Watershed),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

(fig11 <- ggplot(dat_out_full, aes(x = NHD_PctImp2011Cat, y = r_med)) +
    geom_point(alpha = 0.6, size = 3, color = "#A5BA92") +
    labs(x = expression(Percent~Impervious~by~Catchment),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

(fig12 <- ggplot(dat_out_full, aes(x = NHD_PctImp2011Ws, y = r_med)) +
    geom_point(alpha = 0.6, size = 3, color = "#A5BA92") +
    labs(x = expression(Percent~Impervious~by~Watershed),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# quick combination of land use figures
(fig_lu <- (fig9 + fig10) / (fig11 + fig12))

# Effect of Canals
# "95 indicates the least probable interference from a structure of a given type"
(fig13 <- ggplot(dat_out_full, aes(x = Canal, y = r_med)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#8B8D8A", color = "#8B8D8A") +
    labs(x = expression(Likelihood~of~Influence~by~Canals),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# Effect of Dams
# "95 indicates the least probable interference from a structure of a given type"
(fig14 <- ggplot(dat_out_full, aes(x = Dam, y = r_med)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#A7907B", color = "#A7907B") +
    labs(x = expression(Likelihood~of~Influence~by~Dams),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# Nutrients - note, both x axes are LOG SCALED
(fig15 <- ggplot(dat_out_full, aes(x = Nitrate, y = r_med)) +
    geom_point(alpha = 0.6, size = 3, color = "#A6A486") +
    scale_x_log10() +
    labs(x = expression(Mean~Nitrate~(mg/L~NO[3]-N)),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

(fig16 <- ggplot(dat_out_full, aes(x = Orthophosphate, y = r_med)) +
    geom_point(alpha = 0.6, size = 3, color = "#A5BA92") +
    scale_x_log10() +
    labs(x = expression(Mean~OrthoPhosphate~(mg/L~PO[4]-P)),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# Combine figures above.
(fig_r_med <- fig1 + fig1.1 + fig2 +
    fig3.1 + fig9 + fig11 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 2))

# And export for use in the Rmarkdown file.
# ggsave(fig_r_med,
#        filename = "figures/teton_fall22/rmax_6panel_nov.jpg",
#        width = 30,
#        height = 20,
#        units = "cm") # n = 159

(fig_r_supp <- fig3 + fig3.2 + fig7 +
    fig4 + fig4.1 + fig5 +
    fig14 + fig15 + fig16 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 3))

# And export for use in the Rmarkdown file.
# ggsave(fig_r_supp,
#        filename = "figures/teton_fall22/rmax_9panel_nov.jpg",
#        width = 30,
#        height = 30,
#        units = "cm") # n = 159

# Raw GPP and Q for Potomac River site to add alongside CVq figure for job
# application materials:

# Pull out only data of interest
potomac <- dat_in$nwis_01608500

# And filter for year of interest.
potomac12 <- potomac %>%
  filter(date > "2011-12-31") %>%
  filter(date < "2013-01-01")

# And plot using similar format as found in Step 1.
(fig_gpp <- ggplot(potomac12, aes(date, GPP)) +
  geom_point(color="#6CA184", size=3) +
  geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), 
                width=0.2, color="#6CA184") +
  labs(y=expression(atop('GPP', '(g'*~O[2]~m^-2~d^-1*')'))) +
  annotate("rect", xmin = as.Date("2012-02-01"), xmax = as.Date("2012-04-01"),
           ymin = -2, ymax = 12.5, 
           color = "#233D3F", fill = NA, size = 2) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        text = element_text(family = "serif", size = 40)))

(fig_q <- ggplot(potomac12, aes(date, Q)) +
  geom_line(color="#3793EC", size=2) +
  labs(y=expression('Q ('*m^3~s^-1*')'),
       x = "Date") +
  scale_x_date(date_labels = "%b") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(family = "serif", size = 40)))

# Combine figures above.
(fig_potomac <- fig_gpp / fig_q)

# Compile full figure
(figure_app <- ((fig_gpp / fig_q) | fig2.2) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(4,3)))

# ggsave(figure_app,
#        filename = "figures/teton_fall22/gpp_q_r_cvQ.jpg",
#        width = 55,
#        height = 25,
#        units = "cm")

# Plot light availability vs. CVq colored by rmax and mean GPP.
# To compare with Bernhardt et al., 2022 dark & stormy, 
# bright & stable fig.
(fig_schema1 <-ggplot(dat_out_full, aes(x = log10(summerL), 
                                       y = -log10(cvQ),
                                       color = r_med)) +
    geom_point(size = 3) +
    scale_color_viridis() + 
    labs(x = expression(log(Cumulative~Summer~PAR)),
         y = expression(-log(CV[Q])),
         color = expression(r[max])) +
    theme_bw())

(fig_schema2 <-ggplot(dat_out_full, aes(x = log10(summerL), 
                                        y = -log10(cvQ),
                                        color = meanGPP)) +
    geom_point(size = 3) +
    scale_color_viridis(option = "magma") + 
    labs(x = expression(log(Cumulative~Summer~PAR)),
         y = expression(-log(CV[Q])),
         color = expression(GPP)) +
    theme_bw())

(fig_schema <- fig_schema1 + fig_schema2)

# ggsave(fig_schema,
#        filename = "figures/teton_fall22/SummerLight_CVq_r_GPP.jpg",
#        width = 25,
#        height = 10,
#        units = "cm")

#### Value filter for c ####

# Negative c values are not biologically reasonable, so I've 
# removed them.

# First, calculate MEDIAN c values at all the sites that remain
# Using filtered rmax dataset above.
my_159_site_list <- dat_out_full$site_name

dat_out_cmed <- dat_out_df %>%
  filter(site_name %in% my_159_site_list) %>%
  group_by(site_name) %>%
  summarize(c_med = median(c)) %>%
  ungroup()

# And remove negative values.
dat_out_cmed_pos <- dat_out_cmed %>%
  filter(c_med > 0) # Removes 0 sites. Yay!

#### Rhat filter for c ####

# Before proceeding with the analyses, I will be filtering out sites at which
# the model did not converge well for the c parameter.
# Sites with Rhat > 1.05 will not pass muster.

dat_diag_cfilter1 <- dat_diag %>%
  filter(site_name %in% my_159_site_list) %>%
  filter(parameter == "c") %>%
  filter(Rhat < 1.05) # An additional 18 sites drop off.

#### Correllations between s and c values ####

# As an additional model diagnostic, I will evaluate the correlations between
# the s and c values for each iteration at each site and plot them.

# First, I need to filter the original dataset with all iterations by
# the df created above.
my_141_site_list <- dat_diag_cfilter1$site_name

data_out_141 <- dat_out_df %>%
  filter(site_name %in% my_141_site_list)

dat_out_141 <- split(data_out_141, data_out_141$site_name) # remade as list

# Exporting a plot of s vs. c for all iterations for all sites to
# include in the shiny app.

plotting_sc <- function(x) {
  
  names <- unique(x$site_name)
  
  for (i in names){
    
    # create a dataframe at each site
    df <- x %>%
      filter(site_name == i)
    
    # fit a linear model for all iterations
    #fit <- lm(c ~ s, data = df)
    
    # create a plot with r and k for all iterations
    p <- ggplot(df, aes(x = s, y = c)) +
      geom_point(alpha = 0.8) +
      #geom_smooth(method = lm, color = "#E4B3E2") +
      labs(x = "Sensitivity of Persistence Curve (s)",
           y = "Critical Disturbance Threshold (c)") +
           #title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
           #               " p =",signif(summary(fit)$coef[2,4], 5))) +
      theme_bw() +
      theme(text = element_text(size=20))
    
    # display plots
    print(p) 
    
    # save and export plots
    file.name <- paste0("figures/teton_fall22/site_sc_plots/",
                        df$site_name[1],"sc.jpg",sep = "") # create file name
    
    # set specifications for size and resolution of your figure
    ggsave(p,
           filename = file.name,
           width = 8,
           height = 8)
    
  } # close out for loop
  
} # close out function

# test to be sure the function works at a single site
plotting_sc(data_out_141 %>% filter(site_name == "nwis_01124000"))

# And now apply this to the entire dataset.
plotting_sc(data_out_141)

# Roughly, 11 sites had R^2 > 0.5. Lots appear to have log scale relationships.
# Removed lm and model fit print out bc it appeared a bit messy on th app.

#### c Figures ####

# Next, append the positive c values to the Rhat filter to remove
# appropriate sites.
dat_out_cmed_Rhat <- inner_join(dat_diag_cfilter1, dat_out_cmed_pos) 
# 141 sites remaining

# Not removing any sites based on s vs. c plots for the time being.
# REVISIT THIS following co-authors meeting in Nov.

# Now, convert normalized c values to typical discharge values.
dat_maxQ <- dat_in_df %>%
  group_by(site_name) %>%
  summarize(maxQ = max(Q, na.rm = TRUE)) %>%
  ungroup()

dat_together <- left_join(dat_out_cmed_Rhat, dat_maxQ)

dat_together$Qc <- dat_together$c_med*dat_together$maxQ

# And add in 2yr flood to determine Qc:Q2yrf ratio value.

dat_all_together <- left_join(dat_together, dat_2yr)

dat_all_together$Qc_Q2yr <- dat_all_together$Qc/dat_all_together$RI_2yr_Q_cms

# Finally, use dat_in_cvq_L dataset created above for light and CVq.
# And, append this to the larger dataset.
dat_out_yas2 <- left_join(dat_all_together, dat_in_cvq_L)

# Also use dat_site_info dataset created above for site characteristics.
# And append.
dat_out_full_141_1 <- left_join(dat_out_yas2, dat_site_info,
                          by = c("site_name" = "SiteID"))

dat_out_full_141_2 <- left_join(dat_out_full_141_1, dat_site)
dat_out_full_141_3 <- left_join(dat_out_full_141_2, med_width)
dat_out_full_141 <- left_join(dat_out_full_141_3, dat_nuts_w)

# Distribution of Qc/Q2 values:
(fig1qcq2 <- ggplot(dat_out_full_141, aes(x = Qc_Q2yr)) +
    geom_histogram(bins = 60, alpha = 0.8, 
                   fill = "#7E8C69", color = "#7E8C69") +
    labs(x = expression(Q[c]:Q[2~yr]),
         y = "Count") +
    theme_bw())

# CV of Discharge vs. c: note, x axis on LOG SCALE
(fig2qcq2 <- ggplot(dat_out_full_141, aes(x = cvQ, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.8, size = 3,
               color = "#E5A06E") +
    labs(x = expression(CV[Q]),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Mean Daily Light Availability vs. c:
(fig3qcq2 <- ggplot(dat_out_full_141, aes(x = meanL, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.8, size = 3,
               color = "#E6A45A") +
    labs(x = expression(Mean~Daily~PAR),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Stream Order vs. c: Removing singular site w/o order info for now.
(fig4qcq2 <- ggplot(dat_out_full_141 %>%
                  filter(!is.na(Order)), aes(x = Order, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_boxplot(alpha = 0.6, color = "#7E8C69", fill = "#7E8C69") +
    labs(x = expression(Stream~Order),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Stream Width vs. c: note, x axis LOG SCALED
(fig4.1qcq2 <- ggplot(dat_out_full_141, aes(x = width_med, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.6, size = 3, color = "#E7A655") +
    scale_x_log10() + 
    labs(x = expression(Stream~Width~(m)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Latitude vs. c:
(fig5qcq2 <- ggplot(dat_out_full_141, aes(x = Lat_WGS84, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.6, size = 3, color = "#E4957C") +
    labs(x = expression(Latitude),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Longitude vs. c:
(fig6qcq2 <- ggplot(dat_out_full_141, aes(x = Lon_WGS84, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.6, size = 3, color = "#E38678") +
    labs(x = expression(Longitude),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Catchment size vs. c: note, missing Miss. R.
# x-axis also LOG SCALED
(fig7qcq2 <- ggplot(dat_out_full_141, aes(x = NHD_AREASQKM, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.6, size = 3, color = "#E4927B") +
    scale_x_log10() +
    labs(x = expression(Watershed~Area~(km^2)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Land-use vs. c:
(fig8qcq2 <- ggplot(dat_out_full_141, aes(x = LU_category, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_boxplot(alpha = 0.6, color = "#6D4847", fill = "#6D4847") +
    labs(x = expression(Land~Use),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Mean daily GPP vs. c: X axis LOG SCALED
(fig9qcq2 <- ggplot(dat_out_full_141, aes(x = meanGPP, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.8, size = 3,
               color = "#D2A059") +
    scale_x_log10() +
    labs(y = expression(Q[c]:Q[2~yr]),
         x = expression(Mean~Daily~GPP~(gO[2]~m^-2~d^-1))) +
    theme_bw())

# More Land Use vs. c:
(fig10qcq2 <- ggplot(dat_out_full_141, aes(x = NHD_RdDensCat, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.6, size = 3, color = "#CB776D") +
    labs(x = expression(Road~Density~by~Catchment~(km/km^2)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

(fig11qcq2 <- ggplot(dat_out_full_141, aes(x = NHD_RdDensWs, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.6, size = 3, color = "#CB776D") +
    labs(x = expression(Road~Density~by~Watershed),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

(fig12qcq2 <- ggplot(dat_out_full_141, 
                     aes(x = NHD_PctImp2011Cat, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.6, size = 3, color = "#6D4847") +
    labs(x = expression(Percent~Impervious~by~Catchment),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

(fig13qcq2 <- ggplot(dat_out_full_141, 
                     aes(x = NHD_PctImp2011Ws, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.6, size = 3, color = "#6D4847") +
    labs(x = expression(Percent~Impervious~by~Watershed),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Effect of Dams
# "95 indicates the least probable interference from a structure of a given type"
(fig14qcq2 <- ggplot(dat_out_full_141, aes(x = Dam, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#E59D7F", color = "#E59D7F") +
    labs(x = expression(Likelihood~of~Influence~by~Dams),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Combine figures above.
(fig_qcq2 <- fig1qcq2 + fig9qcq2 + fig2qcq2 + 
    fig7qcq2 + fig10qcq2 + fig12qcq2 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 2))

# And export for use in RMarkdown.
# ggsave(fig_qcq2,
#        filename = "figures/teton_fall22/QcQ2_6panel_nov.jpg",
#        width = 30,
#        height = 20,
#        units = "cm") # n = 141

(fig_qcq2_supp <- fig4qcq2 + fig4.1qcq2 + fig14qcq2 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 1))

# ggsave(fig_qcq2_supp,
#        filename = "figures/teton_fall22/QcQ2_3panel_nov.jpg",
#        width = 30,
#        height = 10,
#        units = "cm") # n = 141

#### GPP figures ####

# As of 10/19/22, calculating/exporting RMSE values for the entire dataset
# was proving very computationally intensive. So, for the time being, just
# creating one plot of predicted GPP for inclusion in the coauthors'
# Rmarkdown summary file.

# Going to proceed with site nwis_01649500, NE Anacostia River, because
# it has good data availability and roughly median rmax.
# Re-simulating, because it's simpler than hunting down the right index
# in the matrix above.
Ricker_sim1site <- Ricker_sim_fxn(dat_out$nwis_01649500, dat_in$nwis_01649500)

# And for each day, I would like to calculate
# - mean GPP
# - 97.5% and 2.5% percentiles

# Going to pull out just the predicted GPP values.
data_1site_gpp <- Ricker_sim1site[[1]]

# Calculate median and confidence intervals
median_gpp1 <- apply(data_1site_gpp, 1, median)
lowerci_gpp1 <- apply(data_1site_gpp, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp1 <- apply(data_1site_gpp, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp <- dat_in$nwis_01649500$GPP

# Pull out original dates used
date <- dat_in$nwis_01649500$date

# Bind into a single dataframe
df_1site_pred <- as.data.frame(cbind(median_gpp1, lowerci_gpp1, upperci_gpp1))

df_pred1 <- df_1site_pred %>%
  mutate(date = ymd(date),
         orig_gpp = orig_gpp)

# And finally, calculation the normalized RMSE.
rmse1 <- Ricker_sim1site[[2]]

nRMSE_1site <- nRMSE_fxn(rmse1, dat_in$nwis_01649500)

# And plot
(gpp_plot1 <- ggplot(df_pred1 %>%
                       filter(date > "2010-12-31") %>%
                       filter(date < "2012-01-01"), aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp1), color = "darkolivegreen2", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "NE Anacostia River (MD)") +
    geom_ribbon(aes(ymin = lowerci_gpp1,
                    ymax = upperci_gpp1),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    annotate(geom = "text", x = date("2011-11-01"), y = 10,
             label = paste("nRMSE = ",round(nRMSE_1site, digits = 2)),
             size = 8) + 
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20))) # 2011

# ggsave(gpp_plot1,
#        filename = "figures/teton_fall22/GPP_predGPP_nwis01649500.jpg",
#        width = 40,
#        height = 10,
#        units = "cm")

# End of script.
