## Plotting GPP vs. pred_GPP
## April 25, 2022
## Heili Lowman

# The following script will pull in the initial results of the random effects modeling
# to generate code to plot modeled vs. predicted GPP.

# (the larger model is currently running on Teton, but I want to test things out
# with a like-structured dataset)

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here"), 
       require, character.only=T)

## Source data - with 10 sites
data_in <- readRDS("data_working/stan_10rivers_input_list.rds")
## Model output
data_out <- readRDS("data_working/stan_10rivers_output_Ricker_re_2022_04_22.rds")

# Extract only predicted GPP data from the model
data_out_gpp <- extract(data_out, c("pred_GPP"))
# so, this is an array in the format (rows, columns, matrices)
# in this case, (7500 rows, 192 columns, 10 matrices)

data_out_gpp <- data_out_gpp$pred_GPP

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
#        width = 20,
#        height = 5,
#        units = "cm"
# )

# End of script.
