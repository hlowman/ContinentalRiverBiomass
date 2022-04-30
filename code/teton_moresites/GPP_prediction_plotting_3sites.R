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
data_out <- readRDS("data_working/teton_207rivers_predGPP_withP_021322.rds")
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
data_out_gpp1 <- data_out_gpp[,,1]

# Calculate median and confidence intervals
median_gpp <- apply(data_out_gpp1, 2, median)
lowerci_gpp <- apply(data_out_gpp1, 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upperci_gpp <- apply(data_out_gpp1, 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Pull out original GPP values used
orig_gpp <- data_in$GPP[,1]

# Dummy date column
date <- seq(from = 1, to = 192, by = 1)

# Bind into a single dataframe
df_pred <- as.data.frame(cbind(date, orig_gpp, median_gpp, lowerci_gpp, upperci_gpp))
# And trimming missing dates
df_pred_ed <- df_pred[1:101,]

# And plot
(gpp_plot1 <- ggplot(df_pred_ed, aes(date, orig_gpp)) +
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
          axis.title.y = element_text(size=12)))

# ggsave(("figures/teton_jasm/gpp_vs_pred_gpp_test.png"),
#        width = 12,
#        height = 8,
#        units = "cm"
# )

# End of script.
