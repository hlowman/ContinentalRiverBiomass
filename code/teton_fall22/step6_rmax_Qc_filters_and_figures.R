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
         "calecopal", "viridis", "plotly", "ggbreak"), require, character.only=T)

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

# And the dataset with all HUC delineations.
site_HUC <- readRDS("data_working/HUC12_159sites_120922.rds")

#### Data Prep ####

# Take list containing all input data and make into a df.
dat_in_df <- map_df(dat_in, ~as.data.frame(.x), .id="site_name")

# Take list containing all iterations of parameters and make into a df.
dat_out_df <- map_df(dat_out, ~as.data.frame(.x), .id="site_name")

# Trim HUC dataset.
site_HUC2 <- site_HUC %>%
  dplyr::select(site_name, huc2_id)

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

#### Additional data calculations ####

# Next, append the positive rmax values to the Rhat filter to remove
# appropriate sites.
dat_out_rmed_Rhat <- inner_join(dat_diag_rfilter1, dat_out_rmed_pos) 
# 159 sites remaining

# Not removing any sites based on RMSE values.

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

# Also would like to add site characteristics to this dataset for plotting and linear modeling purposes.
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

# Export for future use.
#saveRDS(dat_out_full, "data_working/rmax_filtered_159sites_113022.rds")

#### rmax Figures ####

# Adding column where the minimum confidence interval of rmax
# truncates at zero.
dat_out_full <- dat_out_full %>%
  mutate(minCI = case_when(`2.5%` < 0 ~ 0,
                           `2.5%` >= 0 ~ `2.5%`))

# Distribution of rmax values:
(fig1 <- ggplot(dat_out_full, aes(x = r_med)) +
  geom_histogram(bins = 60, alpha = 0.8, 
                 fill = "#0B4221", color = "#0B4221") +
  labs(x = expression(Maximum~Growth~Rate~(r[max])),
       y = "Count") +
  theme_bw())

# Mean daily GPP vs. rmax: X axis LOG SCALED
(fig1.1 <- ggplot(dat_out_full, aes(x = meanGPP, y = r_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#285B5D") +
    geom_linerange(alpha = 0.8, 
                   color = "#285B5D",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    scale_x_log10() +
    labs(y = expression(r[max]),
         x = expression(Mean~Daily~GPP~(g~O[2]~m^-2~d^-1))) +
    theme_bw())

# Also, investigating mean daily GPP vs. dams
(fig1.dam <- ggplot(dat_out_full, aes(x = Dam, y = meanGPP)) +
    geom_jitter(alpha = 0.8, size = 3, width = 0.1,
               color = "#285B50") +
    labs(x = "Likelihood of Influence by Dams (%)",
         y = expression(Mean~Daily~GPP~(g~O[2]~m^-2~d^-1))) +
    theme_bw())

# CV of Discharge vs. rmax:
(fig2 <- ggplot(dat_out_full, aes(x = cvQ, y = r_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#518393") +
    geom_linerange(alpha = 0.8, 
                   color = "#518393",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    scale_x_log10() +
    labs(x = expression(CV[Q]),
         y = expression(r[max])) +
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
               color = "#7DACBC") +
    geom_linerange(alpha = 0.8, 
                   color = "#7DACBC",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    scale_x_log10() +
    labs(x = expression(Cumulative~Summer~PAR~(mol~m^-2~d^-1)),
         y = expression(r[max])) +
    theme_bw())

(fig3.2 <- ggplot(dat_out_full, aes(x = summerT, y = r_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#96B0BB") +
    geom_linerange(alpha = 0.8, 
                   color = "#96B0BB",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    labs(x = expression(Mean~Daily~Summer~Temperature~(Celsius)),
         y = expression(r[max])) +
    theme_bw())

# Stream Order vs. rmax: Removing singular site w/o order info for now.
(fig4 <- ggplot(dat_out_full %>%
                  filter(!is.na(Order)), aes(x = Order, y = r_med)) +
    geom_boxplot(alpha = 0.8, color = "black", 
                 fill = "#9BB1BB") +
    labs(x = expression(Stream~Order),
         y = expression(r[max])) +
    theme_bw())

# Stream Width vs. rmax: note, x axis LOG SCALED
(fig4.1 <- ggplot(dat_out_full, aes(x = width_med, y = r_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#C9B7A5") +
    geom_linerange(alpha = 0.8, 
                   color = "#C9B7A5",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    scale_x_log10() + 
    labs(x = expression(Stream~Width~(m)),
         y = expression(r[max])) +
    theme_bw())

# Also, investigating width vs. dams
(fig4.dam <- ggplot(dat_out_full, aes(x = Dam, y = width_med)) +
    geom_jitter(alpha = 0.9, size = 3, width = 0.1,
                color = "#C9B7D9") +
    scale_y_log10() +
    labs(x = "Likelihood of Influence by Dams (%)",
         y = expression(Stream~Width~(m))) +
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
    geom_point(alpha = 0.8, size = 3, color = "#23341E") +
    geom_linerange(alpha = 0.8, 
                   color = "#23341E",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    labs(x = expression(Latitude),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# Longitude vs. rmax:
(fig6 <- ggplot(dat_out_full, aes(x = Lon_WGS84, y = r_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#E4DECE") +
    geom_linerange(alpha = 0.8, 
                   color = "#E4DECE",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    labs(x = expression(Longitude),
         y = expression(r[max])) +
    theme_bw())

# Also did a quick gut check of longtide vs. cvQ
plot(dat_out_full$Lon_WGS84, dat_out_full$cvQ)

# Catchment size vs. rmax: note, missing Miss. R. and x axis LOG SCALED
(fig7 <- ggplot(dat_out_full, aes(x = NHD_AREASQKM, y = r_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#ECBD95") +
    geom_linerange(alpha = 0.8, 
                   color = "#ECBD95",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    scale_x_log10() +
    labs(x = expression(Watershed~Area~(km^2)),
         y = expression(r[max])) +
    theme_bw())

# Land use vs. rmax:
# (fig8 <- ggplot(dat_out_full, aes(x = LU_category, y = r_mean)) +
#     geom_boxplot(alpha = 0.6, color = "#A5BA92", fill = "#A5BA92") +
#     labs(x = expression(Land~Use),
#          y = expression(Maximum~Growth~Rate~(r[max]))) +
#     theme_bw())

(fig9 <- ggplot(dat_out_full, aes(x = NHD_RdDensCat, y = r_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#E9C6A5") +
    geom_linerange(alpha = 0.8, 
                   color = "#E9C6A5",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    labs(x = expression(Road~Density~by~Catchment~(km/km^2)),
         y = expression(r[max])) +
    theme_bw())

(fig10 <- ggplot(dat_out_full, aes(x = NHD_RdDensWs, y = r_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#E9C6A5") +
    geom_linerange(alpha = 0.8, 
                   color = "#E9C6A5",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    labs(x = expression(Road~Density~by~Watershed~(km/km^2)),
         y = expression(r[max])) +
    theme_bw())

(fig11 <- ggplot(dat_out_full, aes(x = NHD_PctImp2011Cat, y = r_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#E4DECE") +
    geom_linerange(alpha = 0.8, 
                   color = "#E4DECE",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    labs(x = expression(Impervious~Land~Cover~by~Catchment~(`%`)),
         y = expression(r[max])) +
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
                 fill = "#A6987F", color = "#A6987F") +
    labs(x = expression(Likelihood~of~Influence~by~Canals),
         y = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# Effect of Dams
# "95 indicates the least probable interference from a structure of a given type"
(fig14 <- ggplot(dat_out_full %>%
                   na.omit(Dam) %>%
                   mutate(DamReOrder = factor(case_when(Dam == "0" ~ "100",
                                                 Dam == "50" ~ "50",
                                                 Dam == "80" ~ "20",
                                                 Dam == "95" ~ "5"),
                                              levels = c("5", "20", "50", "100"))), 
                          aes(x = DamReOrder, y = r_med)) +
    geom_boxplot(alpha = 0.8, 
                 fill = "#E4DECE", color = "black") +
    labs(x = expression(Likelihood~of~Influence~by~Dams~(`%`)),
         y = expression(r[max])) +
    theme_bw())

# An additional figure as I investigate the "dams" category.
(fig14.2 <- ggplot(dat_out_full, aes(x = Dam, y = r_med)) +
    # geom_jitter(alpha = 0.8,
    #             fill = "#0B4229") +
    geom_linerange(alpha = 0.8, position = position_jitter(width = 0.15),
                   color = "#0B4229",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    labs(x = expression(Likelihood~of~Influence~by~Dams~(`%`)),
         y = expression(r[max])) +
    theme_bw())

# Nutrients - note, both x axes are LOG SCALED
(fig15 <- ggplot(dat_out_full, aes(x = Nitrate, y = r_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#0B4221") +
    geom_linerange(alpha = 0.8, 
                   color = "#0B4221",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    scale_x_log10() +
    labs(x = expression(Mean~Nitrate~(mg/L~NO[3]-N)),
         y = expression(r[max])) +
    theme_bw())

(fig16 <- ggplot(dat_out_full, aes(x = Orthophosphate, y = r_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#346575") +
    geom_linerange(alpha = 0.8, 
                   color = "#346575",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    scale_x_log10() +
    labs(x = expression(Mean~OrthoPhosphate~(mg/L~PO[4]-P)),
         y = expression(r[max])) +
    theme_bw())

# Precip
# Annual average precip for local catchment/watershed from HydroATLAS (mm)
site_precip <- site_info %>%
  dplyr::select(SiteID, pre_mm_cyr, pre_mm_uyr)

df_precip_rmax <- left_join(dat_out_full, site_precip, by = c("site_name" = "SiteID"))

(fig17 <- ggplot(df_precip_rmax, aes(x = pre_mm_cyr, y = r_med)) +
    geom_point(alpha = 0.6, size = 3, color = "#79ACBD") +
    geom_linerange(alpha = 0.8, 
                   color = "#79ACBD",
                   aes(ymin = minCI, ymax = `97.5%`)) +
    labs(x = expression(Mean~Annual~Precipitation~by~Catchment~(mm)),
         y = expression(r[max])) +
    theme_bw())

# Combine figures above.
(fig_r_med <- fig1 + fig1.1 + fig2 + fig3.1 +
    fig3.2 + fig4.1 + fig10 + fig14 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 2))

# And export for use in the Rmarkdown file.
# ggsave(fig_r_med,
#        filename = "figures/teton_fall22/rmax_8panel_121522.jpg",
#        width = 40,
#        height = 20,
#        units = "cm") # n = 159

(fig_r_supp <- fig15 + fig16 + fig17 +
    fig4 + fig7 + fig6 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 2))

# And export for use in the Rmarkdown file.
# ggsave(fig_r_supp,
#        filename = "figures/teton_fall22/rmax_6panel_121522.jpg",
#        width = 30,
#        height = 20,
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

# Not removing any sites based on s vs. c plots.

# Now, convert normalized c values to typical discharge values.
dat_maxQ <- dat_in_df %>%
  group_by(site_name) %>%
  summarize(maxQ = max(Q, na.rm = TRUE)) %>%
  ungroup()

dat_together <- left_join(dat_out_cmed_Rhat, dat_maxQ)

# convert both c and confidence interval values.
dat_together$Qc <- dat_together$c_med*dat_together$maxQ
dat_together$Qc2.5 <- dat_together$`2.5%`*dat_together$maxQ
dat_together$Qc97.5 <- dat_together$`97.5%`*dat_together$maxQ

# And add in 2yr flood to determine Qc:Q2yrf ratio value.

dat_all_together <- left_join(dat_together, dat_2yr)

dat_all_together$Qc_Q2yr <- dat_all_together$Qc/dat_all_together$RI_2yr_Q_cms
dat_all_together$Qc_Q2yr2.5 <- dat_all_together$Qc2.5/dat_all_together$RI_2yr_Q_cms
dat_all_together$Qc_Q2yr97.5 <- dat_all_together$Qc97.5/dat_all_together$RI_2yr_Q_cms

# Finally, use dat_in_cvq_L dataset created above for light and CVq.
# And, append this to the larger dataset.
dat_out_yas2 <- left_join(dat_all_together, dat_in_cvq_L)

# Also use dat_site_info dataset created above for site characteristics.
# And append.
dat_out_full_141_1 <- left_join(dat_out_yas2, dat_site_info,
                          by = c("site_name" = "SiteID"))

dat_out_full_141_2 <- left_join(dat_out_full_141_1, dat_site)
dat_out_full_141_3 <- left_join(dat_out_full_141_2, med_width)
dat_out_full_141_4 <- left_join(dat_out_full_141_3, site_HUC2)
dat_out_full_141 <- left_join(dat_out_full_141_4, dat_nuts_w)

dat_out_full_141 <- dat_out_full_141 %>%
# Also creating a new categorical dam column to model by.
mutate(Dam_binary = factor(case_when(
  Dam %in% c("50", "80", "95") ~ "0", # Potential
  Dam == "0" ~ "1", # Certain
  TRUE ~ NA)))

# Export for future use.
#saveRDS(dat_out_full_141, "data_working/QcQ2_filtered_141sites_022823.rds")

# Distribution of c values:
(fig0c <- ggplot(dat_out_full_141, aes(x = c_med)) +
    geom_histogram(bins = 60, alpha = 0.8, 
                   fill = "#262E43", color = "#262E43") +
    labs(x = expression(Critical~Disturbance~Threshold~(Q[c])),
         y = "Count") +
    theme_bw())

# Distribution of Qc/Q2 values:
(fig1qcq2 <- ggplot(dat_out_full_141, aes(x = Qc_Q2yr)) +
    geom_histogram(bins = 60, alpha = 0.8, 
                   fill = "#262E43", color = "#262E43") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    labs(x = expression(Q[c]:Q[2~yr]),
         y = "Count") +
    theme_bw())

# CV of Discharge vs. Qc:Q2: note, x axis on LOG SCALE
(fig2qcq2 <- ggplot(dat_out_full_141, aes(x = cvQ, y = Qc_Q2yr)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#405F8A") +
    geom_linerange(alpha = 0.8, 
                   color = "#405F8A",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    scale_y_log10() + 
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
    geom_boxplot(alpha = 0.6, color = "black", fill = "#5792CC") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() +
    labs(x = expression(Stream~Order),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Stream Width vs. c: note, x axis LOG SCALED
(fig4.1qcq2 <- ggplot(dat_out_full_141, aes(x = width_med, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#1E2F46") +
    geom_linerange(alpha = 0.8, 
                   color = "#1E2F46",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_log10() + 
    scale_y_log10() + 
    labs(x = expression(River~Width~(m)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Latitude vs. c:
(fig5qcq2 <- ggplot(dat_out_full_141, aes(x = Lat_WGS84, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#7E8C69") +
    geom_linerange(alpha = 0.8, 
                   color = "#7E8C69",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(x = expression(Latitude),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Longitude vs. c:
(fig6qcq2 <- ggplot(dat_out_full_141, aes(x = Lon_WGS84, y = Qc_Q2yr)) +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.6, size = 3, color = "#E38678") +
    geom_linerange(alpha = 0.8, 
                   color = "#E38678",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    labs(x = expression(Longitude),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Additional exploratory plot:
(figextra <- ggplot(dat_out_full_141, aes(x = Lon_WGS84, y = cvQ)) +
    geom_point(alpha = 0.6, size = 3, color = "black") +
    labs(x = expression(Longitude),
         y = expression(CV[Q])) +
    theme_bw())

# Catchment size vs. c: note, missing Miss. R.
# x-axis also LOG SCALED
(fig7qcq2 <- ggplot(dat_out_full_141, aes(x = NHD_AREASQKM, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#3B7D6E") +
    geom_linerange(alpha = 0.8, 
                   color = "#3B7D6E",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    scale_y_log10() +
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

# Mean daily GPP vs. Qc:Q2: X axis LOG SCALED
(fig9qcq2 <- ggplot(dat_out_full_141, aes(x = meanGPP, y = Qc_Q2yr)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#486999") +
    geom_linerange(alpha = 0.8, 
                   color = "#486999",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    scale_y_log10() + 
    labs(y = expression(Q[c]:Q[2~yr]),
         x = expression(Mean~Daily~GPP~(gO[2]~m^-2~d^-1))) +
    theme_bw())

# More Land Use vs. c:
(fig10qcq2 <- ggplot(dat_out_full_141, aes(x = NHD_RdDensCat, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#9EB45F") +
    geom_linerange(alpha = 0.8, 
                   color = "#9EB45F",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() + 
    labs(x = expression(Road~Density~by~Catchment~(km/km^2)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

(fig11qcq2 <- ggplot(dat_out_full_141, aes(x = NHD_RdDensWs, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#304969") +
    geom_linerange(alpha = 0.8, 
                   color = "#304969",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() + 
    labs(x = expression(Road~Density~by~Watershed~(km/km^2)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

(fig12qcq2 <- ggplot(dat_out_full_141, 
                     aes(x = NHD_PctImp2011Cat, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#E3907B") +
    geom_linerange(alpha = 0.8, 
                   color = "#E3907B",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
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
(fig14qcq2 <- ggplot(dat_out_full_141 %>%
                       drop_na(Dam_binary), 
                     aes(x = Dam_binary, y = Qc_Q2yr)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#273C57", color = "black") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() + 
    scale_x_discrete(labels = c("5-50%", "100%")) +
    labs(x = expression(Likelihood~of~Influence~by~Dams),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Exploring uncertainty in Qc estimates by dam categories
(fig14qcq2.2 <- ggplot(dat_out_full_141, aes(x = Dam, y = Qc_Q2yr)) +
    geom_linerange(alpha = 0.8, position = position_jitter(width = 0.15),
                   color = "#4D5B90",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() + 
    labs(x = expression(Likelihood~of~Influence~by~Dams),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Effect of Canals
# "95 indicates the least probable interference from a structure of a given type"
(fig15qcq2 <- ggplot(dat_out_full_141, aes(x = Canal, y = Qc_Q2yr)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#E59D7F", color = "#E59D7F") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(x = expression(Likelihood~of~Influence~by~Canals),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Precip
df_precip_Qc_141 <- left_join(dat_out_full_141, site_precip, by = c("site_name" = "SiteID"))

(fig16qcq2 <- ggplot(df_precip_Qc_141, aes(x = pre_mm_cyr, y = Qc_Q2yr)) +
    geom_point(alpha = 0.6, size = 3, color = "#4D5B75") +
    geom_linerange(alpha = 0.8, 
                   color = "#4D5B75",
                   aes(ymin = Qc_Q2yr2.5, ymax = Qc_Q2yr97.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() +
    labs(x = expression(Mean~Annual~Precipitation~by~Catchment~(mm)),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# HUC vs. Qc:Q2: 
(fig17qcq2 <- ggplot(dat_out_full_141, aes(x = huc2_id, y = Qc_Q2yr)) +
    geom_boxplot(alpha = 0.6, color = "black", fill = "#38557A") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() +
    labs(x = expression(Regional~Hydrological~Unit~Code),
         y = expression(Q[c]:Q[2~yr])) +
    theme_bw())

# Combine figures above and export for supplemental figure.
(fig_qcq2_supp <- fig9qcq2 + fig2qcq2 + fig17qcq2 +
    fig11qcq2 + fig14qcq2 + fig4.1qcq2 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 2))

# ggsave(fig_qcq2_supp,
#        filename = "figures/teton_fall22/QcQ2_6panel_022823.jpg",
#        width = 30,
#        height = 20,
#        units = "cm") # n = 141

#### Dam figure ####

# Combining a few of the figures above to examine the effect of/presence
# of Dams.
(fig_dam <- fig14.2 + fig14qcq2.2 +
   fig1.dam + fig4.dam +
   plot_annotation(tag_levels = 'A') +
   plot_layout(nrow = 2))

# ggsave(fig_dam,
#        filename = "figures/teton_fall22/Dams_4panel_031623.jpg",
#        width = 25,
#        height = 20,
#        units = "cm")

#### GPP and NRMSE figures ####

# NOTE - NOT USING nRMSE FILTER FOR FILTERING DATA

# However, I will be creating a figure to demonstrate the universality
# of this modeling approach.

# First, per Phil's suggestion, I will find sites that represent different
# axes of small/large (stream order) and calm/disturbed (CVQ).

(viz_fig <- ggplot(dat_out_full_141, aes(x = Order, y = cvQ)) +
  geom_point(aes(color = c_med, text = site_name)) +
  labs(x = "Stream Order",
       y = "Coefficient of Variation in Discharge") +
  theme_bw())

(viz_plotly <- ggplotly(viz_fig))

(viz_fig2 <- ggplot(dat_out_full, aes(x = Order, y = cvQ)) +
    geom_point(aes(color = r_med, text = site_name)) +
    labs(x = "Stream Order",
         y = "Coefficient of Variation in Discharge") +
    theme_bw())

(viz_plotly2 <- ggplotly(viz_fig2))

(viz_fig3 <- ggplot(dat_out_full_141, aes(x = Order, y = cvQ)) +
    geom_jitter(alpha = 0.8, aes(color = meanGPP, text = site_name)) +
    labs(x = "Stream Order",
         y = "Coefficient of Variation in Discharge") +
    scale_color_viridis() +
    theme_bw())

(viz_plotly3 <- ggplotly(viz_fig3))

(viz_fig4 <- ggplot(dat_out_full_141, aes(x = Order, y = cvQ)) +
    geom_jitter(alpha = 0.8, aes(color = meanL, text = site_name)) +
    labs(x = "Stream Order",
         y = "Coefficient of Variation in Discharge") +
    scale_color_viridis() +
    theme_bw())

(viz_plotly4 <- ggplotly(viz_fig4))

# Create some categories to choose from with axes spanning small/large,
# calm/stormy, and light/dark.

small <- c("1", "2", "3", "4")
big <- c("5", "6", "7", "8")

dat_out_full_141 <- dat_out_full_141 %>%
  mutate(my_groups = factor(case_when(Order %in% small & cvQ < 1 & meanL > 200 ~ "Small_Calm_Light",
                               Order %in% small & cvQ < 1 & meanL < 200 ~ "Small_Calm_Dark",
                               Order %in% small & cvQ > 1 & meanL > 200 ~ "Small_Stormy_Light",
                               Order %in% small & cvQ > 1 & meanL < 200 ~ "Small_Stormy_Dark",
                               Order %in% big & cvQ < 1 & meanL > 200 ~ "Large_Calm_Light",
                               Order %in% big & cvQ < 1 & meanL < 200 ~ "Large_Calm_Dark",
                               Order %in% big & cvQ > 1 & meanL > 200 ~ "Large_Stormy_Light",
                               Order %in% big & cvQ > 1 & meanL < 200 ~ "Large_Stormy_Dark"),
                            levels = c("Small_Calm_Light",
                                       "Small_Calm_Dark",
                                       "Small_Stormy_Light",
                                       "Small_Stormy_Dark",
                                       "Large_Calm_Light",
                                       "Large_Calm_Dark",
                                       "Large_Stormy_Light",
                                       "Large_Stormy_Dark")))

(viz_fig5 <- ggplot(dat_out_full_141 %>% drop_na(Order), aes(x = Order, y = cvQ)) +
    geom_jitter(alpha = 0.8, size = 3, 
                aes(color = my_groups, text = site_name)) +
    labs(x = "Stream Order",
         y = "Coefficient of Variation in Discharge",
         color = "Groups") +
    theme_bw()) ## ooh, this looks good.

# ggsave(viz_fig5,
#        filename = "figures/teton_fall22/Sites_8Groups_111822.jpg",
#        width = 15,
#        height = 10,
#        units = "cm") # n = 141

(viz_plotly5 <- ggplotly(viz_fig5))

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
# Applying function to full 182 site dataset:
# Please note, this takes HOURS to run, so start this early in the day, and
# come back to it.
# Ricker_sim_182sites <- mapply(Ricker_sim_fxn, dat_out, dat_in)

# For some reason, it's yielding two values, so let's see what's happening here.
# Ricker_sim_1site <- Ricker_sim_fxn(dat_out$nwis_01124000, dat_in$nwis_01124000)
# Ok, so it's yielding pred_GPP in the first part and rmse in the second.

# Now, making a longer list to see how it spits out multiple sites to decipher
# my larger output structure.
#dat_in2 <- dat_in[1:2]
#dat_out2 <- dat_out[1:2]

#Ricker_sim_2site <- mapply(Ricker_sim_fxn,dat_out2, dat_in2)
# Viewing this yields nothing, because it's a matrix >_<

# So, for reference:
#predGPP <- Ricker_sim_1site[[1]]
#rmse <- Ricker_sim_1site[[2]]

# But when this is made larger, the portions of the matrix can be accessed by
# indexing by odd and even indices.
# ODD = predGPP
# EVEN = rmse
# So, making a list of even numbers to pull out rmse values.
#my_values <- seq(from = 2, to = 364, by = 2)
#rmse_182sites <- Ricker_sim_182sites[my_values]

# Adding the nRMSE calculation into the function above didn't play nicely with
# the list that existed, so calculating outside instead.
nRMSE_fxn <- function(df, df_orig){
  
  # Calculate the mean RMSE value for each site.
  nRMSE <- mean(df)/(max(df_orig$GPP) - min(df_orig$GPP))
  
}

#nRMSE_182sites <- mapply(nRMSE_fxn, rmse_182sites, dat_in)

# nRMSE_182sitesdf <- as.data.frame(nRMSE_182sites) %>%
#   mutate("site_name" = names(dat_in)) %>%
#   rename("nRMSE" = "nRMSE_182sites")

# Export both sets of results.
# trying to export this first file caused the server to freeze, so only
# exported the nRMSE file for now.
#saveRDS(Ricker_sim_182sites, "data_working/Sim_Ricker_182sites_101922.rds")
#saveRDS(nRMSE_182sitesdf, "data_working/nRMSE_182sites_101922.rds")

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

##### 8-panel NRMSE plot #####

# As of 11/21/22, calculating/exporting RMSE values for 8 sites.

my8sites <- c("nwis_0166818623", "nwis_02217643",
              "nwis_06893350", "nwis_07075250",
              "nwis_05082500", "nwis_13013650",
              "nwis_04176500", "nwis_08374550")

# Trimming input and output datasets for the sites of interest.
dat_out8df <- dat_out_df %>%
  filter(site_name %in% my8sites)

dat_out8 <- split(dat_out8df, dat_out8df$site_name)

dat_in8df <- dat_in_df %>%
  filter(site_name %in% my8sites)

dat_in8 <- split(dat_in8df, dat_in8df$site_name)

# Re-simulating, because it's simpler than hunting down the right index
# in the matrix above. Started ~10:00, Ended ~10:20
Ricker_sim8sites <- mapply(Ricker_sim_fxn, dat_out8, dat_in8)

# And for each day, I would like to calculate
# - median GPP
# - 97.5% and 2.5% percentiles

# Going to pull out just the predicted GPP values.
# So, making a list of odd numbers to pull out predGPP values (see above for reasoning).
my_values <- seq(from = 1, to = 16, by = 2)
data_8site_gpp <- Ricker_sim8sites[my_values]

# Calculate median and confidence intervals
quantile25 <- function(x){quantile(x, probs = 0.025, na.rm = TRUE)}
quantile975 <- function(x){quantile(x, probs = 0.975, na.rm = TRUE)}

pred_gpp8 <- lapply(data_8site_gpp, 
                      function(x) cbind(apply(x, 1, median),
                                        apply(x, 1, quantile25),
                                        apply(x, 1, quantile975)))

# Pull out original GPP values used and sequence #s (for plotting)
orig_gpp_date8 <- lapply(dat_in8, function(x) x %>% select(date, GPP, seq))

# Add names to confidence interval lists
my_names <- c("nwis_0166818623", "nwis_02217643", "nwis_04176500", "nwis_05082500", "nwis_06893350", "nwis_07075250", "nwis_08374550", "nwis_13013650")

names(pred_gpp8) <- my_names
pred_gpp8 <- lapply(pred_gpp8, function(x) as.data.frame(x) %>% 
                      rename("Median" = "V1",
                      "q2.5" = "V2",
                      "q97.5" = "V3")) # OMG YAY!!!!

# Bind into a single dataframe
keys <- unique(c(names(orig_gpp_date8), names(pred_gpp8)))
df_pred8 <- setNames(Map(cbind, orig_gpp_date8[keys], pred_gpp8[keys]), keys)

# And finally, calculate the normalized RMSE.
my_values2 <- seq(from = 2, to = 16, by = 2)
rmse8 <- Ricker_sim8sites[my_values2]

nRMSE_8site <- mapply(nRMSE_fxn, rmse8, dat_in8)

# And plot
# Mill Creek, VA
(gpp_plot8.1 <- ggplot(df_pred8$nwis_0166818623, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Small Order - Steady Flow - High Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2016-06-01"), y = 12.5,
             label = paste("nRMSE = ",round(nRMSE_8site[1], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# Parks Creek, GA
(gpp_plot8.2 <- ggplot(df_pred8$nwis_02217643, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Small Order - Steady Flow - Low Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2016-11-01"), y = 6.5,
             label = paste("nRMSE = ",round(nRMSE_8site[2], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# Tomahawk Creek, KS
(gpp_plot8.3 <- ggplot(df_pred8$nwis_06893350, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Small Order - Turbulent Flow - High Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2012-01-01"), y = 4.5,
             label = paste("nRMSE = ",round(nRMSE_8site[5], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# S. Fork Little Red River, AR
(gpp_plot8.4 <- ggplot(df_pred8$nwis_07075250, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Small Order - Turbulent Flow - Low Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2013-11-01"), y = 11.5,
             label = paste("nRMSE = ",round(nRMSE_8site[6], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# Red River, ND
(gpp_plot8.5 <- ggplot(df_pred8$nwis_05082500, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Large Order - Steady Flow - High Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2016-10-15"), y = 9,
             label = paste("nRMSE = ",round(nRMSE_8site[4], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# Snake River, WY
(gpp_plot8.6 <- ggplot(df_pred8$nwis_13013650, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median, group = seq), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Large Order - Steady Flow - Low Light") +
    scale_x_date(date_labels = "%b %Y") +
    #scale_x_break(c(as.Date("2011-01-01"), as.Date("2014-01-01"))) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, group = seq),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2015-01-01"), y = 30,
             label = paste("nRMSE = ",round(nRMSE_8site[8], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# Raisin River, MI 
(gpp_plot8.7 <- ggplot(df_pred8$nwis_04176500, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Large Order - Turbulent Flow - High Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2014-03-01"), y = 45,
             label = paste("nRMSE = ",round(nRMSE_8site[3], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# Rio Grande, TX
(gpp_plot8.8 <- ggplot(df_pred8$nwis_08374550, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median, group = seq), 
              color = "#609048", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Large Order - Turbulent Flow - Low Light") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, group = seq),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2013-12-01"), y = 7,
             label = paste("nRMSE = ",round(nRMSE_8site[7], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# Combine and export.
(fig_nRMSE <- gpp_plot8.1 + gpp_plot8.2 + 
    gpp_plot8.3 + gpp_plot8.4 +
    gpp_plot8.5 + gpp_plot8.6 +
    gpp_plot8.7 + gpp_plot8.8 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 4))

# ggsave(fig_nRMSE,
#        filename = "figures/teton_fall22/nRMSE_8panel_021023.jpg",
#        width = 40,
#        height = 50,
#        units = "cm") # empty dates fixed :)

##### 5-panel NRMSE plot #####

# As of 02/10/23, calculating/exporting RMSE values for 5
# additional sites that Jud requested.
# Note, need to run data_in_df/data_out_df and 
# Ricker/NRMSE functions above for the below to work.

my5sites <- c("nwis_05524500", "nwis_05515500",
              "nwis_05451210", "nwis_01645704",
              "nwis_0165389205")

# Trimming input and output datasets for the sites of interest.
dat_out5df <- dat_out_df %>%
  filter(site_name %in% my5sites)

dat_out5 <- split(dat_out5df, dat_out5df$site_name)

dat_in5df <- dat_in_df %>%
  filter(site_name %in% my5sites)

dat_in5 <- split(dat_in5df, dat_in5df$site_name)

# Re-simulating, simpler than hunting down the right index
# in the matrix above. Started at 11:00; finished 11:13.
Ricker_sim5sites <- mapply(Ricker_sim_fxn, 
                           dat_out5, dat_in5)

# And for each day, I would like to calculate
# - median GPP
# - 97.5% and 2.5% percentiles

# Going to pull out just the predicted GPP values.
# So, making a list of odd numbers to pull out predGPP values (see above for reasoning).
my_values5 <- seq(from = 1, to = 10, by = 2)
data_5site_gpp <- Ricker_sim5sites[my_values5]

# Use median and confidence intervals' functions from above
pred_gpp5 <- lapply(data_5site_gpp, 
                    function(x) cbind(apply(x, 1, median),
                                      apply(x, 1, quantile25),
                                      apply(x, 1, quantile975)))

# Pull out original GPP values and sequence numbers used
orig_gpp_date5 <- lapply(dat_in5, function(x) x %>% select(date, GPP, Q, seq))

# Add names to confidence interval lists
# be sure that these are in the correct order given output above
my_names5 <- c("nwis_01645704", "nwis_0165389205",
               "nwis_05451210", "nwis_05515500",
               "nwis_05524500")

names(pred_gpp5) <- my_names5
pred_gpp5 <- lapply(pred_gpp5, function(x) as.data.frame(x) %>% 
                      rename("Median" = "V1",
                             "q2.5" = "V2",
                             "q97.5" = "V3")) # OMG YAY!!!!

# Bind into a single dataframe
keys5 <- unique(c(names(orig_gpp_date5), names(pred_gpp5)))
df_pred5 <- setNames(Map(cbind, orig_gpp_date5[keys5], 
                         pred_gpp5[keys5]), keys5)

# And finally, calculate the normalized RMSE.
my_values52 <- seq(from = 2, to = 10, by = 2)
rmse5 <- Ricker_sim5sites[my_values52]

nRMSE_5site <- mapply(nRMSE_fxn, rmse5, dat_in5)

# And plot
# Iroquois River
(gpp_plot8.9 <- ggplot(df_pred5$nwis_05524500, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", linewidth = 1.2) +
    geom_line(aes(date,Q/10),
              color = "blue", linewidth = 1) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Iroquois River, IN",
         subtitle = "Blue line denotes Q/10 (Qmax = 141 cms)") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2016-09-01"), y = 8,
             label = paste("nRMSE = ",round(nRMSE_5site[5], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# Kankakee River, IN
(gpp_plot8.10 <- ggplot(df_pred5$nwis_05515500, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    geom_line(aes(date,Q/10),
              color = "blue", linewidth = 1) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Kankakee River, IN",
         subtitle = "Blue line denotes Q/10 (Qmax = 51 cms)") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2016-10-01"), y = 3,
             label = paste("nRMSE = ",round(nRMSE_5site[4], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# South Fork Iowa River, IA
(gpp_plot8.11 <- ggplot(df_pred5$nwis_05451210 %>%
                          filter(date < "2009-01-01"), 
                        aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median, group = seq), 
              color = "#609048", size = 1.2) +
    geom_line(aes(date,Q/10),
              color = "blue", linewidth = 1) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "S. Fork Iowa River, IA",
         subtitle = "Blue line denotes Q/10 (Qmax = 156 cms)") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, group = seq),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2008-09-01"), y = 15,
             label = paste("nRMSE = ",round(nRMSE_5site[3], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# Difficult Run, VA
(gpp_plot8.12 <- ggplot(df_pred5$nwis_01645704 %>%
                          filter(date > "2014-01-01"), aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median, group = seq), 
              color = "#609048", size = 1.2) +
    geom_line(aes(date,Q),
              color = "blue", linewidth = 1) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Difficult Run, VA",
         subtitle = "Blue line denotes Q (cms)") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, group = seq),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2016-07-01"), y = 6,
             label = paste("nRMSE = ",round(nRMSE_5site[1], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# Accotink Creek, VA
(gpp_plot8.13 <- ggplot(df_pred5$nwis_0165389205, aes(date, GPP)) +
    geom_point(size = 2, color = "#303018") +
    geom_line(aes(date, Median), 
              color = "#609048", size = 1.2) +
    geom_line(aes(date,Q),
              color = "blue", linewidth = 1) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date",
         title = "Accotink Creek, VA",
         subtitle = "Blue line denotes Q (cms)") +
    scale_x_date(date_labels = "%b %Y") +
    geom_ribbon(aes(ymin = q2.5,
                    ymax = q97.5),
                fill = "#90A860", alpha = 0.3) +
    annotate(geom = "text", x = date("2013-10-01"), y = 6,
             label = paste("nRMSE = ",round(nRMSE_5site[2], 
                                            digits = 2)), size = 8) + 
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          axis.title.x = element_text(size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20)))

# Combine and export.
(fig_nRMSE5 <- gpp_plot8.9 + gpp_plot8.10 + 
    gpp_plot8.11 + gpp_plot8.12 +
    gpp_plot8.13 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 3))

# ggsave(fig_nRMSE5,
#        filename = "figures/teton_fall22/nRMSE_5panel_022423.jpg",
#        width = 50,
#        height = 40,
#        units = "cm")

# End of script.
