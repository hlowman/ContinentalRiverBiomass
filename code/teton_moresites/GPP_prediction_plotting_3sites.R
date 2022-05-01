## Plotting GPP vs. pred_GPP
## April 30, 2022
## Heili Lowman

# The following script will pull in the results of the previous modeling
# effort to generate plots of modeled vs. predicted GPP for the JASM poster.

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data - 
data_in <- readRDS("data_working/df_207sites.rds")
## Tidied Teton output for pred GPP only -
data_out207 <- readRDS("data_working/teton_207rivers_model_predGPP_all_iterations_020622.rds")
## Site info -
data_info <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")
## Tidied parameter data -
data_params <- readRDS("data_working/teton_bothmodels_parameters_means_021322.rds") %>%
  filter(model == "with P") %>%
  filter(r_mean > 0)
## Divergences data - 
data_divs <- readRDS("data_working/teton_207rivers_model_divergences_bothmodels_021322.rds")

# what are roughly the 2.5, 50, and 97.5%iles of rmax values?
low <- quantile(data_params$r_mean, probs = 0.025) # 0.02
med <- quantile(data_params$r_mean, probs = 0.5) # 0.12
high <- quantile(data_params$r_mean, probs = .975) # 0.62

# identify sites where divergences are low and rmax aligns close to the above quantiles
data_divs10 <- data_divs %>%
  filter(model == "with P") %>%
  filter(divergences < 10)

data_best <- left_join(data_divs10, data_params)

# 2.5% - nwis_07061270 - East Fork Black River near Lesterville, MO
# 50% - nwis_07332622 - Bois D'Arc Ck at FM 409 nr Honey Grove, TX
# 97.5% - nwis_03025500 - Allegheny River at Franklin, PA

# Use data_out above for pred_GPP values.

# And for each day at each site, I would like to calculate
# - mean GPP
# - 97.5% and 2.5% percentiles

# Going to pull out just the first site
data_out_gpp1 <- data_out207 %>%
  filter(site_name == "nwis_07061270")

data_out_gpp1 <- data_out_gpp1[,-1] # remove site name column
# also checked to be sure data only populated through column 2319
# to match orig_gpp dimensions

# Calculate median and confidence intervals
median_gpp <- apply(data_out_gpp1, 2, median)
lowerci_gpp <- apply(data_out_gpp1, 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp <- apply(data_out_gpp1, 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp <- data_in$nwis_07061270$GPP

# Pull out original dates used
date <- data_in$nwis_07061270$date

# Bind into a single dataframe
df_pred1 <- as.data.frame(cbind(median_gpp, lowerci_gpp, upperci_gpp))

df_pred1 <- df_pred1 %>%
  drop_na() %>%
  mutate(date = ymd(date),
         orig_gpp = orig_gpp)

# And plot
(gpp_plot1 <- ggplot(df_pred1 %>%
                       filter(date > "2014-12-31"), aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp), color = "darkolivegreen2", size = 1.2) +
    labs(y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date") +
    geom_ribbon(aes(ymin = lowerci_gpp,
                    ymax = upperci_gpp),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12))) # 2015-2016?

# And now to pull out just the second site
data_out_gpp2 <- data_out207 %>%
  filter(site_name == "nwis_07332622")

data_out_gpp2 <- data_out_gpp2[,-1] # remove site name column

# Calculate median and confidence intervals
median_gpp2 <- apply(data_out_gpp2, 2, median)
lowerci_gpp2 <- apply(data_out_gpp2, 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp2 <- apply(data_out_gpp2, 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp2 <- data_in$nwis_07332622$GPP

# Pull out original dates used
date2 <- data_in$nwis_07332622$date

# Bind into a single dataframe
df_pred2 <- as.data.frame(cbind(median_gpp2, lowerci_gpp2, upperci_gpp2))

df_pred2 <- df_pred2 %>%
  drop_na() %>%
  mutate(date = ymd(date2),
         orig_gpp = orig_gpp2)

# And plot
(gpp_plot2 <- ggplot(df_pred2 %>%
                       filter(date > "2015-12-31"), aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp2), color = "darkolivegreen2", size = 1.2) +
    labs(#y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date") +
    geom_ribbon(aes(ymin = lowerci_gpp2,
                    ymax = upperci_gpp2),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_blank())) # 2016-2017?

# And now to pull out just the third site
data_out_gpp3 <- data_out207 %>%
  filter(site_name == "nwis_03025500")

data_out_gpp3 <- data_out_gpp3[,-1] # remove site name column

# Calculate median and confidence intervals
median_gpp3 <- apply(data_out_gpp3, 2, median)
lowerci_gpp3 <- apply(data_out_gpp3, 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp3 <- apply(data_out_gpp3, 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp3 <- data_in$nwis_03025500$GPP

# Pull out original dates used
date3 <- data_in$nwis_03025500$date

# Bind into a single dataframe
df_pred3 <- as.data.frame(cbind(median_gpp3, lowerci_gpp3, upperci_gpp3))

df_pred3 <- df_pred3 %>%
  drop_na() %>%
  mutate(date = ymd(date3),
         orig_gpp = orig_gpp3)

# And plot
(gpp_plot3 <- ggplot(df_pred3 %>%
                       filter(date > "2015-12-31"), aes(date, orig_gpp)) +
    geom_point(size = 2, color = "chartreuse4") +
    geom_line(aes(date, median_gpp3), color = "darkolivegreen2", size = 1.2) +
    labs(#y = expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
         x = "Date") +
    geom_ribbon(aes(ymin = lowerci_gpp3,
                    ymax = upperci_gpp3),
                fill = "darkolivegreen2",
                alpha = 0.3) +
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_text(size=12), 
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_blank())) # 2016-2017

gpp_pred_fig <- gpp_plot1 | gpp_plot2 | gpp_plot3

# ggsave(("figures/teton_moresites/gpp_vs_pred_gpp_3site.png"),
#        width = 40,
#        height = 10,
#        units = "cm"
# )

# End of script.
