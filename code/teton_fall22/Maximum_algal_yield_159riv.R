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

# Take list containing all iterations of parameters and make into a df.
dat_out_df <- map_df(dat_out, ~as.data.frame(.x), .id="site_name")

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
            yield_med3 = median(yield_3)) %>%
  ungroup()

# Quick plots.
hist(dat_out_yield_med$yield_med1)
hist(dat_out_yield_med$yield_med2) # still some negative values hmmm...
hist(dat_out_yield_med$yield_med3) # seems similar to formula 1 but x100

# Combine with rmax dataset.
# Note, this will trim down from 182 to 159 sites due to model diagnostics.
dat_yield_rmax <- left_join(dat_rmax, dat_out_yield_med)

# Export for future use.
saveRDS(dat_yield_rmax, "data_working/maxalgalyield_159sites_021323.rds")

#### Figures ####

# Distribution of MAY values:
(fig1 <- ggplot(dat_may_rmax, aes(x = may_med)) +
   geom_histogram(bins = 60, 
                  fill = "#D3E3CA", color = "#D3E3CA") +
   labs(x = expression(Maximum~Algal~Yield),
        y = "Count") +
   theme_bw())

# MAY vs. rmax:
(fig2 <- ggplot(dat_may_rmax, aes(x = r_med, y = may_med)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#CCDFC3") +
    scale_y_log10() +
    labs(y = expression(Maximum~Algal~Yield),
         x = expression(Maximum~Growth~Rate~(r[max]))) +
    theme_bw())

# MAY vs. GPP:
(fig3 <- ggplot(dat_may_rmax, aes(x = meanGPP, y = may_med)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#C7DBBC") +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = expression(Maximum~Algal~Yield),
         x = expression(Mean~Daily~GPP~(gO[2]~m^-2~d^-1))) +
    theme_bw())

# MAY vs. cvQ:
(fig4 <- ggplot(dat_may_rmax, aes(x = cvQ, y = may_med)) +
    geom_point(alpha = 0.9, size = 3,
               color = "#C1D7B6") +
    scale_y_log10() +
    labs(x = expression(CV[Q]),
         y = expression(Maximum~Algal~Yield)) +
    theme_bw())

# MAY vs. summer light:
(fig5 <- ggplot(dat_may_rmax, aes(x = summerL, y = may_med)) +
    geom_point(alpha = 0.8, size = 3,
               color = "#B7CFAC") +
    scale_y_log10() +
    labs(x = expression(Cumulative~Summer~PAR~(mol~m^-2~d^-1)),
         y = expression(Maximum~Algal~Yield)) +
    theme_bw())

# MAY vs. stream width:
(fig6 <- ggplot(dat_may_rmax, aes(x = width_med, y = may_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#ABC1A0") +
    scale_x_log10() + 
    scale_y_log10() +
    labs(x = expression(Stream~Width~(m)),
         y = expression(Maximum~Algal~Yield)) +
    theme_bw())

# MAY vs. longitude:
(fig7 <- ggplot(dat_may_rmax, aes(x = Lon_WGS84, y = may_med)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#9EB393") +
    scale_y_log10() +
    labs(x = expression(Longitude),
         y = expression(Maximum~Algal~Yield)) +
    theme_bw())

# MAY vs. Road density:
(fig8 <- ggplot(dat_may_rmax, aes(x = NHD_RdDensCat, y = may_med)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#92A587") +
    scale_y_log10() +
    labs(x = expression(Road~Density~by~Catchment~(km/km^2)),
         y = expression(Maximum~Algal~Yield)) +
    theme_bw())

# MAY vs. NO3:
(fig9 <- ggplot(dat_may_rmax, aes(x = Nitrate, y = may_med)) +
    geom_point(alpha = 0.8, size = 3, color = "#7D8D70") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Mean~Nitrate~(mg/L~NO[3]-N)),
         y = expression(Maximum~Algal~Yield)) +
    theme_bw())

# MAY vs. PO4:
(fig10 <- ggplot(dat_may_rmax, aes(x = Orthophosphate, y = r_med)) +
    geom_point(alpha = 0.8, size = 3, 
               color = "#687659") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression(Mean~OrthoPhosphate~(mg/L~PO[4]-P)),
         y = expression(Maximum~Algal~Yield)) +
    theme_bw())

# MAY vs. dams:
(fig11 <- ggplot(dat_may_rmax, aes(x = Dam, y = may_med)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#545F43", color = "#545F43") +
    scale_y_log10() +
    labs(x = expression(Likelihood~of~Influence~by~Dams),
         y = expression(Maximum~Algal~Yield)) +
    theme_bw())

# MAY vs. canals:
(fig12 <- ggplot(dat_may_rmax, aes(x = Canal, y = may_med)) +
    geom_boxplot(alpha = 0.6, 
                 fill = "#464F35", color = "#464F35") +
    scale_y_log10() +
    labs(x = expression(Likelihood~of~Influence~by~Canals),
         y = expression(Maximum~Algal~Yield)) +
    theme_bw())

# Combine figures above.
(fig_may_med_log <- fig1 + fig2 + fig3 +
    fig4 + fig5 + fig6 +
    fig7 + fig8 + fig9 +
    fig10 + fig11 + fig12 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(nrow = 4))

# And export for use in the Rmarkdown file.
# ggsave(fig_may_med_log,
#        filename = "figures/teton_fall22/maxalgyield_12panel_120222.jpg",
#        width = 30,
#        height = 40,
#        units = "cm") # n = 159

# End of script.
