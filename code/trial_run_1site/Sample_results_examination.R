# The following code was used in the original iteration of the Biomass1_FitModels.R
# script before repurposing/editing it for use on Teton.

# I've copied the code below so that I have a record of how to extract
# parameters for later analysis.

## View
# launch_shinystan(PM_outputlist_AR$nwis_08447300)
PM_outputlist_Ricker <- readRDS("data_working/stan_1riv_output_Ricker_2021_07_07.rds")
launch_shinystan(PM_outputlist_Ricker$nwis_01608500)

PM_outputlist_Ricker_noP <- readRDS("data_working/stan_1riv_output_RickernoP_2021_07_07.rds")
launch_shinystan(PM_outputlist_Ricker_noP$nwis_01608500)
# launch_shinystan(PM_outputlist_Ricker_tvr$nwis_11044000)

# Extract parameters
mparams <- extract(PM_outputlist_Ricker$nwis_01608500, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o"))

mparams_noP <- extract(PM_outputlist_Ricker_noP$nwis_01608500, c("r","lambda","B","pred_GPP","sig_p","sig_o"))

# Take the median of biomass, P, and predicted GPP outputs
biomass_output <- as.data.frame(mparams$B) %>%
  pivot_longer(cols = everything(), names_to = "day") %>%
  group_by(day) %>%
  summarize(median_b = median(value)) %>%
  ungroup()

persistence_output <- as.data.frame(mparams$P) %>%
  pivot_longer(cols = everything(), names_to = "day") %>%
  group_by(day) %>%
  summarize(median_p = median(value)) %>%
  ungroup()

pred_gpp_output <- as.data.frame(mparams$pred_GPP) %>%
  pivot_longer(cols = everything(), names_to = "day") %>%
  group_by(day) %>%
  summarize(median_gpp = median(value)) %>%
  ungroup()

# create dataset of initial data used
rawdat <- data.frame(df$nwis_01608500)

outputs <- left_join(biomass_output, persistence_output) %>%
  left_join(pred_gpp_output) %>%
  # need to make day numeric (if not labeled properly)
  mutate(day = str_remove(day, "[V]")) %>%
  mutate(day = as.numeric(day)) %>%
  # sort by newly created day column
  arrange(day) %>%
  # add calendar dates
  mutate(date_ymd = rawdat$date) %>%
  # and pivot for facetting
  pivot_longer(cols = c(median_b, median_p, median_gpp),
               names_to = "parameter")

# examine structure of outputs
str(outputs)

# Plot daily medians of parameters
outputs_plot <- ggplot(outputs, aes(x = date_ymd, y = value)) +
  geom_point(aes(color = parameter)) +
  scale_color_viridis(discrete = TRUE) +
  facet_grid(rows = vars(parameter), scales = "free") +
  labs(x = "Date",
       y = "Median daily value",
       title = "South Branch Potomac River (WV)",
       subtitle = "Ricker Model Results - 2015") +
  theme_bw() +
  theme(legend.position = "none")

outputs_plot

# ggsave(plot = outputs_plot,
#        filename = "figures/ricker_results/nwis_01608500_2015.jpg",
#        width = 8,
#        height = 6)

# Repeat with model results from removing persistence term.
# Take the median of biomass, P, and predicted GPP outputs
biomass_output_noP <- as.data.frame(mparams_noP$B) %>%
  pivot_longer(cols = everything(), names_to = "day") %>%
  group_by(day) %>%
  summarize(median_b = median(value)) %>%
  ungroup()

pred_gpp_output_noP <- as.data.frame(mparams_noP$pred_GPP) %>%
  pivot_longer(cols = everything(), names_to = "day") %>%
  group_by(day) %>%
  summarize(median_gpp = median(value)) %>%
  ungroup()

outputs_noP <- left_join(biomass_output_noP, pred_gpp_output_noP) %>%
  # need to make day numeric (if not labeled properly)
  mutate(day = str_remove(day, "[V]")) %>%
  mutate(day = as.numeric(day)) %>%
  # sort by newly created day column
  arrange(day) %>%
  # add calendar dates
  mutate(date_ymd = rawdat$date) %>%
  # and pivot for facetting
  pivot_longer(cols = c(median_b, median_gpp),
               names_to = "parameter")

# examine structure of outputs
str(outputs_noP)

# Plot daily medians of parameters
outputs_plot_noP <- ggplot(outputs_noP, aes(x = date_ymd, y = value)) +
  geom_point(aes(color = parameter)) +
  scale_color_viridis(discrete = TRUE) +
  facet_grid(rows = vars(parameter), scales = "free") +
  labs(x = "Date",
       y = "Median daily value",
       title = "South Branch Potomac River (WV)",
       subtitle = "Ricker Model Results Without Persistence Term - 2015") +
  theme_bw() +
  theme(legend.position = "none")

outputs_plot_noP

# ggsave(plot = outputs_plot_noP,
#        filename = "figures/ricker_results/nwis_01608500_noP_2015.jpg",
#        width = 8,
#        height = 6)

# Also per Joanna's suggestion, going to plot results on top of one another to
# see how things match up.

outputs1 <- outputs %>%
  # choose only gpp predictions
  filter(parameter == "median_gpp") %>%
  select(date_ymd, value) %>%
  rename(date = date_ymd, 
         GPP = value) %>%
  # add column to denote model structure
  mutate(model = "With P")

outputs2 <- outputs_noP %>%
  # choose only gpp predictions
  filter(parameter == "median_gpp") %>%
  select(date_ymd, value) %>%
  rename(date = date_ymd, 
         GPP = value) %>%
  # add column to denote model structure
  mutate(model = "Without P")

outputs3 <- rawdat %>%
  select(date, GPP) %>%
  mutate(model = "Raw data")

# bind all datasets together for plotting
alloutputs <- rbind(outputs1, outputs2, outputs3)

# Plot GPP for both models and raw data
(comparison_plot <- ggplot(alloutputs, aes(x = date, y = GPP)) +
    geom_point(aes(color = model)) +
    scale_color_viridis(discrete = TRUE) +
    facet_grid(rows = vars(model), scales = "free") +
    labs(x = "Date",
         y = "GPP",
         color = " ",
         title = "South Branch Potomac River (WV)",
         subtitle = "Comparison of Appling Data and Ricker Model Results - 2015") +
    theme_bw()+
    theme(legend.position = "none"))

# ggsave(plot = comparison_plot,
#        filename = "figures/ricker_results/comparison_nwis_01608500_2015.jpg",
#        width = 8,
#        height = 6)
