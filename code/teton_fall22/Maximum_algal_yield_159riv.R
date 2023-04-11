## Estimating maximum algal yields
## December 1, 2022
## Heili Lowman

# The following script will calculate the maximum algal yields.
# Then, it will examine if they have any relationship with other covariates.

# Edits on February 13th, 2023 : I've added in two additional metrics of yield
# to try out, so the three formulations are as follows.

# (1) Yield = ((-0.5*rmax/lambda)*(e^0.5*rmax)) + (0.5*rmax/lambda) [based
# on calculations provided by C. Yakulic]
# (2) Yield = (rmax*(0.5-(0.07*rmax))) / -lambda [based on Scheuerell 2016 
# PeerJ]
# (3) Yield = (lambert_W0(exp(1-rmax))) / -lambda [using the 'gsl' package]

#### Setup ####

# Load packages.
lapply(c("lubridate","tidyverse", "here", "viridis",
         "reshape2","ggExtra","patchwork", "gsl"), require, character.only=T)

# Load datasets.

# Load in list containing all iterations of site-level parameters.
dat_out <- readRDS("data_working/teton_182rivers_model_params_all_iterations_101522.rds")

# And remaining covariate data for plotting.
dat_rmax <- readRDS("data_working/rmax_filtered_159sites_113022.rds")

# And the dataset with all HUC delineations.
site_HUC <- readRDS("data_working/HUC12_159sites_120922.rds")

# And the dataset with site metadata.
site <- fread("data_raw/site_data.tsv")

# And the updated dataset with nutrient data.
nuts <- readRDS("data_working/USGS_WQP_nuts_aggsite_022322.rds")

# And the dataset containing exceedances.
exc <- readRDS("data_working/Qc_exceedances_159sites_021423.rds")

# Take list containing all iterations of parameters and make into a df.
dat_out_df <- map_df(dat_out, ~as.data.frame(.x), .id="site_name")

# Trim HUC dataset.
site_HUC2 <- site_HUC %>%
  dplyr::select(site_name, huc2_id)

# Trim nutrient dataset only for P.
nuts_P <- nuts %>%
  filter(CharacteristicName == "Phosphorus") %>%
  rename(Phosphorus = "mean_mg_L")

#### Yield Formulas ####

# Equation to calculate max. algal yield using the first formula is:
# yield = (-0.5*r*exp(0.5*r)/lambda) + (0.5*r/lambda)

# Equation to calculate max. algal yield using the second formula is:
# yield = (r*(0.5-(0.07*r)))/-lambda

# Equation to calculate max algal yield using the third formula is:
# yield = (lambert_W0(exp(1-r)))/-lambda

# Create new columns to calculate yield for each iteration at each site.
dat_out_df <- dat_out_df %>%
  mutate(yield_1 = (-0.5*r*exp(0.5*r)/lambda) + (0.5*r/lambda)) %>%
  mutate(yield_2 = (r*(0.5-(0.07*r)))/-lambda) %>%
  mutate(yield_3 = (lambert_W0(exp(1-r)))/-lambda)

# Note, I've added an extra negative value to the denominator in the second
# and third formulas since our lambda values are negative, but the manuscript
# assumes b values are positive (see Table 1, Scheuerell 2016 PeerJ).

# Calculate median yield values for each site and each yield metric.
dat_out_yield_med <- dat_out_df %>%
  group_by(site_name) %>%
  summarize(yield_med1 = median(yield_1),
            yield_med2 = median(yield_2),
            yield_2.5_2 = quantile(yield_2, probs = 0.025),
            yield_97.5_2 = quantile(yield_2, probs = 0.975),
            yield_med3 = median(yield_3),
            lambda_med = median(lambda)) %>%
  ungroup() %>%
  # Some lower quantiles of yield were negative, so need to
  # filter those to become 0 instead so they plot.
  mutate(yield_2.5_2ed = case_when(yield_2.5_2 < 0 ~ 0,
                                   TRUE ~ yield_2.5_2))

# Quick plots.
hist(dat_out_yield_med$yield_med1)
hist(dat_out_yield_med$yield_med2) # still some negative values hmmm...
# Oh, this is because we don't yet have the poorly performing sites
# filtered out.
hist(dat_out_yield_med$yield_med3) # seems similar to formula 1 but x100

# For ease, making the original list into a df.
dat_in_df <- map_df(dat_in, ~as.data.frame(.x), .id="site_name")

# And count years in each record at each site.
dat_years <- dat_in_df %>%
  group_by(site_name) %>%
  summarize(years = n_distinct(year)) %>%
  ungroup()

# Combine with rmax dataset.
# Note, this will trim down from 182 to 159 sites due to model diagnostics.
dat_yield_rmax_1 <- left_join(dat_rmax, dat_out_yield_med)
# Combine with HUC dataset
dat_yield_rmax_2 <- left_join(dat_yield_rmax_1, site_HUC2)
# Edit nutrient dataset so it has a matching sitename column
nuts_P$site_name <- str_replace_all(nuts_P$MonitoringLocationIdentifier, 'USGS-', 'nwis_')
# And combine primary df with nutrient data
dat_yield_rmax_3 <- left_join(dat_yield_rmax_2, nuts_P)
# And combine with exceedances data.
dat_yield_rmax_4 <- left_join(dat_yield_rmax_3, exc)
# and combine with years data to calculate exceedances/year.
dat_yield_rmax_5 <- left_join(dat_yield_rmax_4, dat_years)

dat_yield_rmax <- dat_yield_rmax_5 %>%
  # Create new column for annual exceedance values used in models.
  mutate(exc_ev_y = total_exc_events/years) %>%
  # Also creating a new categorical dam column to model by.
  mutate(Dam_binary = factor(case_when(
    Dam %in% c("50", "80", "95") ~ "0", # Potential
    Dam == "0" ~ "1", # Certain
    TRUE ~ NA)))

# Export for future use.
# saveRDS(dat_yield_rmax, "data_working/maxalgalyield_159sites_041123.rds")

#### Figures ####

# Going to make three large plots of these formulae to see how 
# they compare in terms of covariate trends also.

##### Formula #1: #####
# Distribution of MAY values:
(fig1.1 <- ggplot(dat_yield_rmax, aes(x = yield_med1)) +
   geom_histogram(bins = 60, 
                  fill = "#D3E3CA", color = "#C7DBBC") +
   labs(x = expression(Maximum~Algal~Yield~(Yakulic/Lowman)),
        y = "Count") +
   theme_bw())

# MAY vs. rmax:
(fig1.2 <- ggplot(dat_yield_rmax, aes(x = r_med, y = yield_med1)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#B7CFAC") +
    scale_y_log10() +
    labs(y = expression(Maximum~Algal~Yield~(Yakulic/Lowman)),
         x = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# MAY vs. lambda:
(fig1.2.2 <- ggplot(dat_yield_rmax, aes(x = lambda_med, y = yield_med1)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#687659") +
    scale_y_log10() +
    labs(y = expression(Maximum~Algal~Yield~(Yakulic/Lowman)),
         x = expression(Lambda)) +
    theme_bw())

# MAY vs. GPP:
(fig1.3 <- ggplot(dat_yield_rmax, aes(x = meanGPP, y = yield_med1)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#C7DBBC") +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = expression(Maximum~Algal~Yield~(Yakulic/Lowman)),
         x = expression(Mean~Daily~GPP~(gO[2]~m^-2~d^-1))) +
    theme_bw())

# MAY vs. cvQ:
(fig1.4 <- ggplot(dat_yield_rmax, aes(x = cvQ, y = yield_med1)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#C1D7B6") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(CV[Q]),
         y = expression(Maximum~Algal~Yield~(Yakulic/Lowman))) +
    theme_bw())

# MAY vs. summer light:
(fig1.5 <- ggplot(dat_yield_rmax, aes(x = summerL, y = yield_med1)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#B7CFAC") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Cumulative~Summer~PAR~(mol~m^-2~d^-1)),
         y = expression(Maximum~Algal~Yield~(Yakulic/Lowman))) +
    theme_bw())

# MAY vs. stream width:
(fig1.6 <- ggplot(dat_yield_rmax, aes(x = width_med, y = yield_med1)) +
    geom_point(alpha = 0.8, size = 3, color = "#ABC1A0") +
    scale_x_log10() + 
    scale_y_log10() +
    labs(x = expression(Stream~Width~(m)),
         y = expression(Maximum~Algal~Yield~(Yakulic/Lowman))) +
    theme_bw())

# MAY vs. longitude:
(fig1.7 <- ggplot(dat_yield_rmax, aes(x = Lon_WGS84, y = yield_med1)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#9EB393") +
    scale_y_log10() +
    labs(x = expression(Longitude),
         y = expression(Maximum~Algal~Yield~(Yakulic/Lowman))) +
    theme_bw())

# MAY vs. Road density:
(fig1.8 <- ggplot(dat_yield_rmax, aes(x = NHD_RdDensCat, y = yield_med1)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#92A587") +
    #scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Road~Density~by~Catchment~(km/km^2)),
         y = expression(Maximum~Algal~Yield~(Yakulic/Lowman))) +
    theme_bw())

# MAY vs. NO3:
(fig1.9 <- ggplot(dat_yield_rmax, aes(x = Nitrate, y = yield_med1)) +
    geom_point(alpha = 0.8, size = 3, color = "#7D8D70") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Mean~Nitrate~(mg/L~NO[3]-N)),
         y = expression(Maximum~Algal~Yield~(Yakulic/Lowman))) +
    theme_bw())

# MAY vs. PO4:
(fig1.10 <- ggplot(dat_yield_rmax, aes(x = Orthophosphate, 
                                       y = yield_med1)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#687659") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Mean~OrthoPhosphate~(mg/L~PO[4]-P)),
         y = expression(Maximum~Algal~Yield~(Yakulic/Lowman))) +
    theme_bw())

# MAY vs. dams:
(fig1.11 <- ggplot(dat_yield_rmax, aes(x = Dam, y = yield_med1)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#545F43", color = "#545F43") +
    scale_y_log10() +
    labs(x = expression(Likelihood~of~Influence~by~Dams),
         y = expression(Maximum~Algal~Yield~(Yakulic/Lowman))) +
    theme_bw())

# MAY vs. canals:
(fig1.12 <- ggplot(dat_yield_rmax, aes(x = Canal, y = yield_med1)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#464F35", color = "#464F35") +
    scale_y_log10() +
    labs(x = expression(Likelihood~of~Influence~by~Canals),
         y = expression(Maximum~Algal~Yield~(Yakulic/Lowman))) +
    theme_bw())

# Combine figures above.
(fig_yield_med1 <- fig1.1 + fig1.2 + fig1.3 +
    fig1.4 + fig1.5 + fig1.6 +
    fig1.7 + fig1.8 + fig1.9 +
    fig1.10 + fig1.11 + fig1.12 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 4))

# And export for use in the Rmarkdown file.
# ggsave(fig_yield_med1,
#        filename = "figures/teton_fall22/maxalgyield1_12panel_021323.jpg",
#        width = 30,
#        height = 40,
#        units = "cm") # n = 159

# And making a smaller, three-panel plot.
(fig_yield_med1_tiny <- fig1.1 + fig1.2 + fig1.2.2 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 1))

# ggsave(fig_yield_med1_tiny,
#        filename = "figures/teton_fall22/maxalgyield1_3panel_021323.jpg",
#        width = 30,
#        height = 10,
#        units = "cm") # n = 159

##### Formula #2: #####
# Distribution of MAY values:
(fig2.1 <- ggplot(dat_yield_rmax, aes(x = yield_med2)) +
    geom_histogram(bins = 60, 
                   fill = "#F6EECF", color = "#EEDCB4") +
    labs(x = expression(Maximum~Algal~Yield~(Scheuerell)),
         y = "Count") +
    theme_bw())

# MAY vs. rmax:
(fig2.2 <- ggplot(dat_yield_rmax, aes(x = r_med, y = yield_med2)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#D0B692") +
    scale_y_log10() +
    labs(y = expression(Maximum~Algal~Yield~(Scheuerell)),
         x = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# MAY vs. lambda:
(fig2.2.2 <- ggplot(dat_yield_rmax, aes(x = lambda_med, y = yield_med2)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#53261B") +
    scale_y_log10() +
    labs(y = expression(Maximum~Algal~Yield~(Scheuerell)),
         x = expression(Lambda)) +
    theme_bw())

# MAY vs. GPP:
(fig2.3 <- ggplot(dat_yield_rmax, aes(x = meanGPP, y = yield_med2)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#96AA8B") +
    geom_linerange(alpha = 0.8, 
                   color = "#96AA8B",
                   aes(ymin = yield_2.5_2ed, ymax = yield_97.5_2)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = expression(a[max]),
         x = expression(Mean~Daily~GPP~(gO[2]~m^-2~d^-1))) +
    theme_bw())

# MAY vs. cvQ:
(fig2.4 <- ggplot(dat_yield_rmax, aes(x = cvQ, y = yield_med2)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#8A9C7E") +
    geom_linerange(alpha = 0.8, 
                   color = "#8A9C7E",
                   aes(ymin = yield_2.5_2ed, ymax = yield_97.5_2)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(CV[Q]),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. summer light:
(fig2.5 <- ggplot(dat_yield_rmax, aes(x = summerL, y = yield_med2)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#7B8B6E") +
    geom_linerange(alpha = 0.8, 
                   color = "#7B8B6E",
                   aes(ymin = yield_2.5_2ed, ymax = yield_97.5_2)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Cumulative~Summer~PAR~(mol~m^-2~d^-1)),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. stream width:
(fig2.6 <- ggplot(dat_yield_rmax, aes(x = width_med, y = yield_med2)) +
    geom_point(alpha = 0.8, size = 3, color = "#454F34") +
    geom_linerange(alpha = 0.8, 
                   color = "#454F34",
                   aes(ymin = yield_2.5_2ed, ymax = yield_97.5_2)) +
    scale_x_log10() + 
    scale_y_log10() +
    labs(x = expression(Stream~Width~(m)),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. longitude:
(fig2.7 <- ggplot(dat_yield_rmax, aes(x = Lon_WGS84, y = yield_med2)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#A27E65") +
    scale_y_log10() +
    labs(x = expression(Longitude),
         y = expression(Maximum~Algal~Yield~(Scheuerell))) +
    theme_bw())

# MAY vs. Road density:
(fig2.8 <- ggplot(dat_yield_rmax, 
                  aes(x = NHD_RdDensWs, y = yield_med2)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#3A422D") +
    geom_linerange(alpha = 0.8, 
                   color = "#3A422D",
                   aes(ymin = yield_2.5_2ed, ymax = yield_97.5_2)) +
    #scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Road~Density~by~Watershed~(km/km^2)),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. NO3:
(fig2.9 <- ggplot(dat_yield_rmax, aes(x = Nitrate, y = yield_med2)) +
    geom_point(alpha = 0.8, size = 3, color = "#343B29") +
    geom_linerange(alpha = 0.8, 
                   color = "#343B29",
                   aes(ymin = yield_2.5_2ed, ymax = yield_97.5_2)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Mean~Nitrate~(mg/L~NO[3]-N)),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. P:
(fig2.10 <- ggplot(dat_yield_rmax, aes(x = Phosphorus, 
                                       y = yield_med2)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#2F3525") +
    geom_linerange(alpha = 0.8, 
                   color = "#2F3525",
                   aes(ymin = yield_2.5_2ed, ymax = yield_97.5_2)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Mean~Dissolved~Phosphorus~(mg/L~P)),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. dams:
(fig2.11 <- ggplot(dat_yield_rmax %>%
                     drop_na(Dam_binary), 
                   aes(x = Dam_binary, y = yield_med2)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#4D583C", color = "black") +
    scale_y_log10() +
    scale_x_discrete(labels = c("5-50%", "100%")) +
    labs(x = expression(Likelihood~of~Influence~by~Dams),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. canals:
(fig2.12 <- ggplot(dat_yield_rmax, aes(x = Canal, y = yield_med2)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#291611", color = "#291611") +
    scale_y_log10() +
    labs(x = expression(Likelihood~of~Influence~by~Canals),
         y = expression(Maximum~Algal~Yield~(Scheuerell))) +
    theme_bw())

# MAY vs. summer temp:
(fig2.13 <- ggplot(dat_yield_rmax, aes(x = summerT, y = yield_med2)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#5C694C") +
    geom_linerange(alpha = 0.8, 
                   color = "#5C694C",
                   aes(ymin = yield_2.5_2ed, ymax = yield_97.5_2)) +
    scale_y_log10() +
    labs(x = expression(Mean~Summer~Temperature~(degree~C)),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. HUC: 
(fig2.14 <- ggplot(dat_yield_rmax, aes(x = huc2_id, y = yield_med2)) +
    geom_boxplot(alpha = 0.6, color = "black", fill = "#6C7A5D") +
    scale_y_log10() +
    labs(x = expression(Regional~Hydrological~Unit~Code),
         y = expression(a[max])) +
    theme_bw())

# MAY vs. annual exceedances:
(fig2.15 <- ggplot(dat_yield_rmax, aes(x = exc_ev_y, y = yield_med2)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#6C7A5D") +
    geom_linerange(alpha = 0.8, 
                   color = "#6C7A5D",
                   aes(ymin = yield_2.5_2ed, ymax = yield_97.5_2)) +
    #scale_x_log10() +
    scale_y_log10() +
    labs(y = expression(a[max]),
         x = expression(Annual~Exceedances~of~Q[c])) +
    theme_bw())

# Combine figures above.
(fig_yield_med2 <- fig2.3 + fig2.4 + fig2.5 + fig2.15 +
    fig2.14 + fig2.13 + fig2.11 + fig2.6 + fig2.8 +
    fig2.9 + fig2.10 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 3))

# And export for use in the Rmarkdown file.
# ggsave(fig_yield_med2,
#        filename = "figures/teton_fall22/maxalgyield2_11panel_041123.jpg",
#        width = 40,
#        height = 30,
#        units = "cm") # n = 159

# And combine and export tiny figure as well.
(fig_yield_med2_tiny <- fig2.1 + fig2.2 + fig2.2.2 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 1))

# ggsave(fig_yield_med2_tiny,
#        filename = "figures/teton_fall22/maxalgyield2_3panel_021323.jpg",
#        width = 30,
#        height = 10,
#        units = "cm") # n = 159

##### Formula #3: #####
# Distribution of MAY values:
(fig3.1 <- ggplot(dat_yield_rmax, aes(x = yield_med3)) +
   geom_histogram(bins = 60, 
                  fill = "#A1CAF6", color = "#75A1DE") +
   labs(x = expression(Maximum~Algal~Yield~(Lambert)),
        y = "Count") +
   theme_bw())

# MAY vs. rmax:
(fig3.2 <- ggplot(dat_yield_rmax, aes(x = r_med, y = yield_med3)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#5982BD") +
    scale_y_log10() +
    labs(y = expression(Maximum~Algal~Yield~(Lambert)),
         x = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# MAY vs. lambda:
(fig3.2.2 <- ggplot(dat_yield_rmax, aes(x = lambda_med, y = yield_med3)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#304969") +
    scale_y_log10() +
    labs(y = expression(Maximum~Algal~Yield~(Lambert)),
         x = expression(Lambda)) +
    theme_bw())

# MAY vs. GPP:
(fig3.3 <- ggplot(dat_yield_rmax, aes(x = meanGPP, y = yield_med3)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#75A1DE") +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = expression(Maximum~Algal~Yield~(Lambert)),
         x = expression(Mean~Daily~GPP~(gO[2]~m^-2~d^-1))) +
    theme_bw())

# MAY vs. cvQ:
(fig3.4 <- ggplot(dat_yield_rmax, aes(x = cvQ, y = yield_med3)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#628ED1") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(CV[Q]),
         y = expression(Maximum~Algal~Yield~(Lambert))) +
    theme_bw())

# MAY vs. summer light:
(fig3.5 <- ggplot(dat_yield_rmax, aes(x = summerL, y = yield_med3)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#5982BD") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Cumulative~Summer~PAR~(mol~m^-2~d^-1)),
         y = expression(Maximum~Algal~Yield~(Lambert))) +
    theme_bw())

# MAY vs. stream width:
(fig3.6 <- ggplot(dat_yield_rmax, aes(x = width_med, y = yield_med3)) +
    geom_point(alpha = 0.8, size = 3, color = "#5075AA") +
    scale_x_log10() + 
    scale_y_log10() +
    labs(x = expression(Stream~Width~(m)),
         y = expression(Maximum~Algal~Yield~(Lambert))) +
    theme_bw())

# MAY vs. longitude:
(fig3.7 <- ggplot(dat_yield_rmax, aes(x = Lon_WGS84, y = yield_med3)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#486999") +
    scale_y_log10() +
    labs(x = expression(Longitude),
         y = expression(Maximum~Algal~Yield~(Lambert))) +
    theme_bw())

# MAY vs. Road density:
(fig3.8 <- ggplot(dat_yield_rmax, aes(x = NHD_RdDensCat, y = yield_med3)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#405F8A") +
    #scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Road~Density~by~Catchment~(km/km^2)),
         y = expression(Maximum~Algal~Yield~(Lambert))) +
    theme_bw())

# MAY vs. NO3:
(fig3.9 <- ggplot(dat_yield_rmax, aes(x = Nitrate, y = yield_med3)) +
    geom_point(alpha = 0.8, size = 3, color = "#38557A") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Mean~Nitrate~(mg/L~NO[3]-N)),
         y = expression(Maximum~Algal~Yield~(Lambert))) +
    theme_bw())

# MAY vs. PO4:
(fig3.10 <- ggplot(dat_yield_rmax, aes(x = Orthophosphate, 
                                       y = yield_med3)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#304969") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Mean~OrthoPhosphate~(mg/L~PO[4]-P)),
         y = expression(Maximum~Algal~Yield~(Lambert))) +
    theme_bw())

# MAY vs. dams:
(fig3.11 <- ggplot(dat_yield_rmax, aes(x = Dam, y = yield_med3)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#273C57", color = "#273C57") +
    scale_y_log10() +
    labs(x = expression(Likelihood~of~Influence~by~Dams),
         y = expression(Maximum~Algal~Yield~(Lambert))) +
    theme_bw())

# MAY vs. canals:
(fig3.12 <- ggplot(dat_yield_rmax, aes(x = Canal, y = yield_med3)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#1E2F46", color = "#1E2F46") +
    scale_y_log10() +
    labs(x = expression(Likelihood~of~Influence~by~Canals),
         y = expression(Maximum~Algal~Yield~(Lambert))) +
    theme_bw())

# Combine figures above.
(fig_yield_med3 <- fig3.1 + fig3.2 + fig3.3 +
    fig3.4 + fig3.5 + fig3.6 +
    fig3.7 + fig3.8 + fig3.9 +
    fig3.10 + fig3.11 + fig3.12 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 4))

# And export for use in the Rmarkdown file.
# ggsave(fig_yield_med3,
#        filename = "figures/teton_fall22/maxalgyield3_12panel_021323.jpg",
#        width = 30,
#        height = 40,
#        units = "cm") # n = 159

# And again create smaller three-paneled plot and export.
(fig_yield_med3_tiny <- fig3.1 + fig3.2 + fig3.2.2 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 1))

# ggsave(fig_yield_med3_tiny,
#        filename = "figures/teton_fall22/maxalgyield3_3panel_021323.jpg",
#        width = 30,
#        height = 10,
#        units = "cm") # n = 159

#### accrual quantiles ####
# 95th %tile and above
quantile(dat_yield_rmax$yield_med2, probs = c(0.95)) # 10.91745
a95above <- dat_yield_rmax %>%
  filter(yield_med2 >= 10.91745) # 8 sites
a95above_sites <- site %>%
  filter(site_name %in% a95above$site_name) # shares 5 of 8 sites with rmax

# 5th %tile and below
quantile(dat_yield_rmax$yield_med2, probs = c(0.05)) # 0.353083
a5below <- dat_yield_rmax %>%
  filter(yield_med2 <= 0.353083) # 8 sites
a5below_sites <- site %>%
  filter(site_name %in% a5below$site_name) # shares 3 of 8 sites with rmax

# Exploring commonalities
# GPP
mean(a95above$meanGPP) # 8.098374
mean(a5below$meanGPP) # 0.3311079

# Size
mean(a95above$width_med) # 51.33848
mean(a5below$width_med) # 7.801826

# Discharge
mean(a95above$cvQ) # 1.100216
mean(a5below$cvQ) # 1.83351

# End of script.
