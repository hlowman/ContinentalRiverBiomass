## 207 rivers data output from Teton
## Created: February 6, 2022
## Heili Lowman

# The following code will do some preliminary examination of the output
# of the 207 site dataset sent to Teton and run through the Ricker model
# earlier this week.

# We'll be loading in data generated having run through the "TetonOutput_
# processing.R" script on the Pinyon server.

# The first part of the code will include only the output of the version
# with the P term included. The output for sites where the P term has
# been removed will be added later.

# Load packages
lapply(c("calecopal", "cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", "here",
         "ggrepel", "patchwork", "grid","gridExtra"), require, character.only=T)

#### Data Import ####

# Load dataset loaded into Teton.
data_in <- readRDS("data_working/df_207sites.rds")

# Load dataset with site information.
data_info <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")

# Load dataset with parameters' data (all iterations).
data_params <- readRDS("data_working/teton_207rivers_model_parameters_all_iterations_020622.rds")
data_params2 <- readRDS("data_working/teton_131rivers_model_parameters_all_iterations_021322.rds")

# Load dataset with model divergences for each site.
data_divs <- readRDS("data_working/teton_207rivers_model_divergences_020622.rds")
data_divs2 <- readRDS("data_working/teton_131rivers_model_divergences_021322.rds")

# Load dataset with model diagnostics for each site.
data_diags <- readRDS("data_working/teton_207rivers_model_diagnostics_020622.rds")
data_diags2 <- readRDS("data_working/teton_131rivers_model_diagnostics_021322.rds")

#### Data Tidying ####

# First, need to join all of the datasets above.

# all iterations of the parameters
data_params <- data_params %>%
  mutate(model = "with P") # add model designation

data_params2 <- data_params2 %>%
  mutate(model = "without P",
         s = NA,
         c = NA) %>% # add model designation and missing columns
  select(site_name, r, lambda, s, c, k, model) # and reorder

data_params_all <- rbind(data_params, data_params2) # and join together

# Export dataset
# saveRDS(data_params_all, file = "data_working/teton_207rivers_model_parameters_all_iterations_bothmodels_021322.rds")

# all divergences
data_divs <- data_divs %>%
  mutate(model = "with P") # add model designation

data_divs2 <- data_divs2 %>%
  mutate(model = "without P") # add model designation

data_divs_all <- rbind(data_divs, data_divs2) # and join together

# Export dataset
# saveRDS(data_divs_all, file = "data_working/teton_207rivers_model_divergences_bothmodels_021322.rds")

# all diagnostics
data_diags <- data_diags %>%
  mutate(model = "with P") # add model designation

data_diags2 <- data_diags2 %>%
  mutate(model = "without P") # add model designation

data_diags_all <- rbind(data_diags, data_diags2) # and join together

# Export dataset
# saveRDS(data_diags_all, file = "data_working/teton_207rivers_model_diagnostics_bothmodels_021322.rds")



# And now to calculate means by site.
data_params_means <- data_params %>%
  group_by(site_name) %>%
  summarize(r_mean = mean(r),
            k_mean = mean(k),
            s_mean = mean(s),
            c_mean = mean(c))

# And now to bind the values with site attributes.
data_together <- left_join(data_params_means, data_info, by = "site_name")

# Export dataset
# saveRDS(data_together, file = "data_working/teton_207rivers_model_parameters_means_020622.rds")

####      Figures         ####

#### Examining Relationship between r & k ####

# Mean r vs. mean k.
# Exporting a plot of r vs. k for all iterations for all sites to
# include in the shiny app (similar to the covariate plots).

plotting_rk <- function(x) {
  
  names <- unique(x$site_name)
  
  for (i in names){
  
  # create a dataframe at each site
  df <- x %>%
    filter(site_name == i)
  
  # join with site information
  #df <- left_join(df.1, data_info, by = "site_name")
  
  # create a plot with r and k for all iterations
  p <- ggplot(df, aes(x = r, y = logK)) +
       geom_point() +
       labs(x = "Maximum Growth Rate (r)",
              y = "Log of Carrying Capacity (K)") +
       theme_bw() +
       theme(text = element_text(size=20))
    
    # display plots
    print(p) 
    
    # save and export plots
    file.name <- paste0("figures/teton_moresites/site_rk_plots/",
                        df$site_name[1],"rk.jpg",sep = "") # create file name
    
    # set specifications for size and resolution of your figure
    ggsave(p,
           filename = file.name,
           width = 8,
           height = 8)
    
  } # close out for loop
  
} # close out function

data_params <- data_params %>%
  mutate(logK = log10(k))# add log(k) column

# test to be sure the function works at a single site
plotting_rk(data_params %>% filter(site_name == "nwis_01124000"))

# And now apply this to the entire dataset.
# plotting_rk(data_params)

#### Max Growth Rate / Carrying Capacity ####

# Basic plot of r values vs. stream order.
fig1 <- data_together %>%
  filter(r_mean > 0) %>%
  #filter(k_mean > 0) %>%
  mutate(so = factor(NHD_STREAMORDE)) %>%
  ggplot(aes(x = so, y = r_mean, fill = so)) +
  scale_fill_manual(values = cal_palette("sbchannel", n = 9, type = "continuous")) + # custom colors
  geom_boxplot(alpha = 0.9) +
  geom_jitter(color = "black", alpha = 0.5, width = 0.1) +
  labs(x = "Stream Order",
       y = "Maximum Growth Rate (r)") +
  theme_bw() +
  theme(text = element_text(size=12), legend.position = "none")

fig1

# ggsave(plot = fig1,
#        filename = "figures/teton_moresites/r_strord.jpg",
#        width = 6,
#        height = 4)

fig2 <- data_together %>%
  #filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  mutate(so = factor(NHD_STREAMORDE)) %>%
  mutate(logK = log10(k_mean)) %>%
  ggplot(aes(x = so, y = logK, fill = so)) +
  scale_fill_manual(values = cal_palette("sbchannel", n = 9, type = "continuous")) +
  geom_boxplot(alpha = 0.9) +
  geom_jitter(color = "black", alpha = 0.5, width = 0.1) +
  labs(x = "NHD Stream Order",
       y = "Log of Carrying Capacity (K)") +
  theme_bw() +
  theme(text = element_text(size=12), legend.position = "none")

fig2

# ggsave(plot = fig2,
#        filename = "figures/teton_moresites/k_strord.jpg",
#        width = 6,
#        height = 4)

# Combine the above boxplots into a single figure.
fig1.2 <- (fig1 | fig2) +
  plot_annotation(tag_levels = "A")

fig1.2

# ggsave(plot = fig1.2,
#        filename = "figures/teton_moresites/r_k_strord.jpg",
#        width = 10,
#        height = 6)

# A quick pairs plot.
fig_pairs <- data_together %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  select(r_mean, k_mean, s_mean, c_mean, NHD_STREAMORDE) %>%
  pairs()

# Plotting overall r and k means for all sites.
fig3a <- data_together %>%
  filter(r_mean > 0) %>% # removes 10 sites
  filter(k_mean > 0) %>% # removes another 4 sites
  mutate(logK = log10(k_mean)) %>%
  ggplot(aes(x = r_mean, y = logK)) +
  geom_point(shape = 21, size = 4, alpha = 0.75) +
  labs(x = "Maximum Growth Rate (r)",
       y = "Log of Carrying Capacity (K)") +
  theme_bw() +
  theme(text = element_text(size=12))

fig3a

# ggsave(plot = fig3a,
#        filename = "figures/teton_moresites/r_K_allsites.jpg",
#        width = 6,
#        height = 6)

# Create another version of this figure colored by stream order
fig3b <- data_together %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  mutate(logK = log10(k_mean)) %>%
  mutate(so = factor(NHD_STREAMORDE)) %>%
  ggplot(aes(x = r_mean, y = logK, 
             fill = so)) +
  geom_point(shape = 21, size = 4, alpha = 0.9) +
  scale_fill_manual(values = cal_palette("sbchannel", n = 9, type = "continuous")) + # custom colors
  labs(x = "Maximum Growth Rate (r)",
       y = "Log of Carrying Capacity (K)") +
  theme_bw() +
  theme(text = element_text(size=12))

fig3b

#### STOPPED HERE ON FEBRUARY 6 ####

#### Critical Discharge / Sensitivity of Persistence Curve ####

# and exploring disturbance metrics
fig4a <- ggplot(data_together, aes(x = c_mean, y = s_mean, 
                                   fill = site_name, label = site_name)) +
  geom_point(shape = 21, size = 4, alpha = 0.75) +
  labs(x = "Critical Discharge (c)",
       y = "Sensitivity of Persistence Curve (s)",
       title = "Full Dataset") +
  geom_text_repel(size=3) +
  theme_bw() +
  theme(legend.position = "none")

fig4a

# Removing negative r and K values.
fig4b <- data_together %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  ggplot(aes(x = c_mean, y = s_mean, 
             fill = site_name, label = site_name)) +
  geom_point(shape = 21, size = 5, alpha = 0.75) +
  labs(x = "Critical Discharge (c)",
       y = "Sensitivity of Persistence Curve (s)") +
  scale_fill_manual(values = cal_palette("creek", n = 28, type = "continuous")) + # custom colors
  geom_text_repel(data = subset(data_together, c_mean < 0.25 & s_mean > 250), size = 4) +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")

fig4b

# Panel by stream order with negative growth parameters removed
fig4c <- data_together %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  ggplot(aes(x = c_mean, y = s_mean, fill = site_name)) +
  geom_point(shape = 21, size = 4, alpha = 0.75) +
  labs(x = "Critical Discharge (c)",
       y = "Steepness of Persistence Curve (s)",
       title = "Paneled by Stream Order - Negative Values Removed") +
  facet_wrap(~NHD_STREAMORDE, scales = "free", ncol = 4)+
  theme_bw() +
  theme(legend.position = "none")

fig4c

# Creating some boxplots for s and c as well

# Basic plot of s values vs. stream order.
fig2.1 <- data_together %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  mutate(so = factor(NHD_STREAMORDE)) %>%
  ggplot(aes(x = so, y = s_mean, fill = so)) +
  scale_fill_manual(values = cal_palette("sbchannel", n = 7, type = "continuous")) + # custom colors
  geom_boxplot(alpha = 0.75) +
  geom_jitter(color = "black", alpha = 0.5, width = 0.1) +
  labs(x = "NHD Stream Order",
       y = "Sensitivity of Persistence Curve (s)") +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")

fig2.1

# Basic plot of c values vs. stream order.
fig2.2 <- data_together %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  mutate(so = factor(NHD_STREAMORDE)) %>%
  ggplot(aes(x = so, y = c_mean, fill = so)) +
  scale_fill_manual(values = cal_palette("sbchannel", n = 7, type = "continuous")) + # custom colors
  geom_boxplot(alpha = 0.75) +
  geom_jitter(color = "black", alpha = 0.5, width = 0.1) +
  labs(x = "NHD Stream Order",
       y = "Critical Discharge (c)") +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")

fig2.2

# creating composite figure for export
fig2_full <- fig2.1 + fig2.2 + fig4b +
  plot_annotation(tag_levels = 'A')

fig2_full

# ggsave(fig2_full,
#        filename = "figures/teton_34sites/s_c_paneled.jpg",
#        width = 15,
#        height = 5)

# Adding figures for Modelscape demonstration

# Site that performed poorly - Reedy Creek, FL nwis_02266300
# Site that performed well - South Branch Potomac, WV nwis_01608500
(figc1 <- data_out_params_df %>%
  filter(site_name == "nwis_02266300" | site_name == "nwis_01608500") %>%
  mutate(site = factor(site_name, levels = c("nwis_02266300", "nwis_01608500"))) %>%
  ggplot(aes(x = c, color = site, fill = site)) +
  geom_histogram(alpha = 0.75) +
  scale_color_manual(values = cal_palette("sierra2")) + # custom colors
  scale_fill_manual(values = cal_palette("sierra2")) + # custom colors
  geom_density(alpha = 0.75) +
  labs(x = "Critical Discharge (c)",
       y = "Density") +
  facet_wrap(.~site_name, scales = "free") +
  theme_bw() +
  theme(text = element_text(size=12), legend.position = "none"))

(figs1 <- data_out_params_df %>%
    filter(site_name == "nwis_02266300" | site_name == "nwis_01608500") %>%
    mutate(site = factor(site_name, levels = c("nwis_02266300", "nwis_01608500"))) %>%
    ggplot(aes(x = s, color = site, fill = site)) +
    geom_histogram(alpha = 0.75) +
    scale_color_manual(values = cal_palette("sierra2")) + # custom colors
    scale_fill_manual(values = cal_palette("sierra2")) + # custom colors
    geom_density(alpha = 0.75) +
    labs(x = "Sensitivity of Persistence (s)",
         y = "Density") +
    facet_wrap(.~site_name, scales = "free") +
    theme_bw() +
    theme(text = element_text(size=12), legend.position = "none"))

# creating composite figure for export
figcs_full <- figc1 / figs1

figcs_full

# ggsave(figcs_full,
#        filename = "figures/presentations/s_c_histograms.jpg",
#        width = 10,
#        height = 10)

#### Persistence Curves ####

# Adapted from Joanna's code in: 
# RiverBiomass/code/Biomass5_Fig_PersistenceCurves.R

## Extract and summarize parameters
par_Ricker <- lapply(data_out, function(x) rstan::extract(x, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o")))

## mean parameter function
par_mean <- function(par) {
  ## Find the mean
  # of all iterations for all parameters
  mean_par <- lapply(par, function(x) mean(x))
  
  # of each iteration for GPP, P, and B parameters
  mean_pred_GPP_ts <- apply(par$pred_GPP,2,mean)
  mean_P_ts <- apply(par$P,2,mean)
  mean_B_ts <- apply(par$B,2,mean)
  
  ## Compile in list and return
  mean_par_ts <- list(mean_par, mean_pred_GPP_ts, mean_P_ts, mean_B_ts)
  names(mean_par_ts) <- c("par","pred_GPP","P","B")
  return(mean_par_ts)
}

# Apply to dataset of extracted parameters
meanpar_R <- lapply(par_Ricker, function(x) par_mean(x))

## Plot persistence
# function to pull out information of interest (i.e., discharge and persistence data)
persistence_list <- function(y, data){
  Ppars <- list()
  
  for(i in 1:length(y)){
    # treat each site as a new dataframe
    df <- y[[i]]
    # same for input data, each site is a new df
    dat <- data[[i]]
    # pull out info re: discharge (Q), critical discharge (c), & sensitivity of the curve (s)
    Ppars[[i]] <- list("tQ"=dat$tQ,"range"=range(dat$tQ),
                       "c"=df$c,"s"=df$s, 
                       "site_name"=dat$site_name[1])
  }
  
  names(Ppars) <- names(y)
  
  return(Ppars)
  
}

# And now map this to the output parameter data and the input raw data
P_R <- persistence_list(par_Ricker, data_in)

## plot
# function calculating 97.5%, 50%, and 2.5% percentiles
# of the persistence term (temp)
plotting_P_dat <- function(x){
  pq <- seq(x$range[1],x$range[2], length=length(x$s))
  p_median <- numeric()
  p_up <- numeric()
  p_down <- numeric()
  name <- substring(deparse(substitute(x)),7)
  for(i in 1:length(pq)){
    temp <- exp(-exp(x$s*(pq[i]-x$c)))
    p_median[i] <- median(temp)
    p_up[i] <- quantile(temp, probs = 0.975)
    p_down[i] <- quantile(temp, probs = 0.025)
  }
  df <- as.data.frame(as.matrix(cbind(pq,p_median, p_up, p_down)))
  df$site_name <- x$site_name
  
  return(df)
}

# use percentile-calculating function on aggregated persistence/discharge dataset
# create in the step above
P_dat_R1 <- lapply(P_R, function(z) plotting_P_dat(z))

#P_dat_R1$PM <- "Ricker" # removed since the Ricker model is the only kind of model I'm running

# rename said file
P_df <- P_dat_R1

Persistence_plots <- function(site, df, site_info, P_df){
  
  Q_sub <- df[[site]]
  #Q_sub$pq <- Q_sub$tQ
  Q_sub$p_for_q <- 0.5
  
  ## critical Q based on velocity # removed for now since I don't have bankfull discharge data
  # at all 34 sites, and already commented out of plot created below
  #crit_Q <- site_info[which(site_info$site_name == site),]$RI_2yr_Q
  
  ## convert relativized Q to original values
  P <- P_df[[site]] # changed this line based on troubleshooting below
  P$Q <- P$pq*max(Q_sub$Q, na.rm = T)
  
  ## critical Q based on GPP - Q correction needed
  c <- meanpar_R[[site]]$par$c*max(Q_sub$Q, na.rm = T)
  
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  ## Plot creation
  # use percentile dataset, discharge on x, persistence values on y
  Persist_plot <- ggplot(P, aes(Q, p_median))+
    scale_x_continuous(trans = "log", labels = scaleFUN)+
    geom_point(data=Q_sub, aes(Q, p_for_q), color="white")+ # ???
    geom_line(size=1.5, alpha=0.9, color="chartreuse4")+ # median p values
    geom_ribbon(data=P, aes(ymin=p_down, ymax=p_up), alpha=0.3, fill="chartreuse4", color=NA)+ # confidence intervals of p values
    theme(panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12, angle=45, hjust=1),
          # axis.title = element_blank(), # want to keep axis labels in for now
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15))+
    labs(x="Range of Standardized Discharge",y="Persistence", title = P$site_name[1])+
    scale_y_continuous(limits=c(0,1))+
    geom_vline(xintercept = c, size=1, linetype="dashed")
  
  
  Persist_plot2 <- ggExtra::ggMarginal(Persist_plot, data=Q_sub, type="histogram",
                                       size=4, x = Q, margins = "x", color="black",
                                       fill="deepskyblue4", xparams = list(alpha=0.8))
  
  # display combined plot
  print(Persist_plot2)
  
  # save and export plots
  file.name <- paste0("figures/teton_34sites/site_persistence_curves/",
                      Q_sub$site_name[1],"_persist.jpg",sep = "") # create file name
  
  # set specifications for size and resolution of your figure
  ggsave(Persist_plot2,
         filename = file.name,
         width = 8,
         height = 8)
  
}

# create the full list of 34 sites
site_list <- levels(as.factor(data_info$site_name))

# apply function above to create plots at each site
plots <- lapply(site_list, function(x) Persistence_plots(x,data_in,data_info,P_df))
# Joanna uses the structure: function(site, df, site_info, P_df)
# which caused me some confusion initially, since it's sourcing data files from other scripts
# in the same repository. Instead, I've chosen to name the files the same as they are named
# in the import steps at the very beginning of *this* script, for consistency.
# Keep this in mind when troubleshooting future scripts!

#### Persistence Curve Troubleshooting ####
# So, this function isn't working, so going to try working through it step by step below, at a single site:
Q_sub <- data_in[["nwis_01608500"]]
Q_sub$p_for_q <- 0.5

## convert relativized Q to original values
#P <- P_df[which(P_df$site_name == "nwis_01608500"),] # here was the error I was getting before
P <- P_df[["nwis_01608500"]] # trying this instead to see if it works the same...it seems to.
P$Q <- P$pq*max(Q_sub$Q, na.rm = T)

## critical Q based on GPP - Q correction needed
# critical discharge * maximum discharge at a site
c <- meanpar_R[["nwis_01608500"]]$par$c*max(Q_sub$Q, na.rm = T)

scaleFUN <- function(x) sprintf("%.1f", x)

## Plot
Persist_plot <- ggplot(P, aes(Q, p_median))+
  scale_x_continuous(trans = "log", labels = scaleFUN)+
  geom_point(data=Q_sub, aes(Q, p_for_q), color="white")+
  geom_line(size=1.5, alpha=0.9, color="chartreuse4")+
  geom_ribbon(data=P, aes(ymin=p_down, ymax=p_up), alpha=0.3, fill="chartreuse4", color=NA)+
  theme(panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle=45, hjust=1),
        #axis.title = element_blank(), 
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  labs(x="Range of Standardized Discharge",y="Persistence", title = "nwis_01608500")+
  scale_y_continuous(limits=c(0,1))+
  geom_vline(xintercept = c, size=1, linetype="dashed")


Persist_plot2 <- ggExtra::ggMarginal(Persist_plot, data=Q_sub, type="histogram",
                                     size=4, x = Q, margins = "x", color="black",
                                     fill="deepskyblue4", xparams = list(alpha=0.8))

Persist_plot
Persist_plot2 # WOOT!


names(plots) <- site_list

#### All Parameter Comparisons ####

# Now, to compare growth and disturbance parameters
# r vs. c
fig5a <- data_together %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  ggplot(aes(x = c_mean, y = r_mean, 
             fill = site_name, label = site_name)) +
  geom_point(shape = 21, size = 5, alpha = 0.75) +
  scale_fill_manual(values = cal_palette("creek", n = 28, type = "continuous")) + # custom colors
  labs(x = "Critical Discharge (c)",
       y = "Maximum Growth Rate (r)") +
  geom_text_repel(data = subset(data_together, r_mean > 0.3), size = 4) +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")

fig5a

# r vs. s
fig5b <- data_together %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  ggplot(aes(x = s_mean, y = r_mean, 
             fill = site_name, label = site_name)) +
  geom_point(shape = 21, size = 5, alpha = 0.75) +
  scale_fill_manual(values = cal_palette("creek", n = 28, type = "continuous")) +
  labs(x = "Sensitivity of Persistence Curve (s)",
       y = "Maximum Growth Rate (r)") +
  geom_text_repel(data = subset(data_together, r_mean > 0.3 | s_mean > 300), size = 4) +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")

fig5b

# k vs. c
fig5c <- data_together %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  ggplot(aes(x = c_mean, y = k_mean, 
             fill = site_name, label = site_name)) +
  geom_point(shape = 21, size = 5, alpha = 0.75) +
  scale_fill_manual(values = cal_palette("creek", n = 28, type = "continuous")) +
  labs(x = "Critical Discharge (c)",
       y = "Carrying Capacity (K)") +
  geom_text_repel(data = subset(data_together, k_mean > 30), size = 4) +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")

fig5c

# k vs. s
fig5d <- data_together %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  ggplot(aes(x = s_mean, y = k_mean, 
             fill = site_name, label = site_name)) +
  geom_point(shape = 21, size = 5, alpha = 0.75) +
  scale_fill_manual(values = cal_palette("creek", n = 28, type = "continuous")) +
  labs(x = "Sensitivity of Persistence Curve (s)",
       y = "Carrying Capacity (K)") +
  geom_text_repel(data = subset(data_together, k_mean > 30 | s_mean > 300), size = 4) +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")

fig5d

full_fig5 <- (fig5a + fig5b) / (fig5c + fig5d) +
  plot_annotation(tag_levels = 'A')

full_fig5

# ggsave(full_fig5,
#        filename = "figures/teton_34sites/rk_vs_cs.jpg",
#        width = 11,
#        height = 11)

#### Light vs. Growth Parameters ####

# Create a function to calculate mean light availability
calc_mean_light <- function(df){
  mean(df$PAR_new) # this column either pulls light from Appling or PAR_surface from Savoy
}

# And now map this to the entire output list.
data_light <- map(data_in, calc_mean_light)

# And create a dataframe
data_light_params <- map_df(data_light, ~as.data.frame(.x), .id="site_name") %>%
  rename(light_mean = `.x`) %>%
  left_join(data_together, by = "site_name")

# Light vs.maximum growth rate
# light vs. r
fig_light1 <- data_light_params %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  mutate(light_cat = factor(light_mean)) %>% # adding column for coloration purposes
  ggplot(aes(x = light_mean, y = r_mean, fill = light_cat)) +
  geom_point(shape = 21, size = 5, alpha = 0.75) +
  scale_fill_manual(values = cal_palette("canary", n = 28, type = "continuous")) + # custom colors
  labs(x = "Light Availability (ppfd?)",
       y = "Maximum Growth Rate (r)") +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")

fig_light1

# Light vs.carrying capacity
# light vs. K
fig_light2 <- data_light_params %>%
  filter(r_mean > 0) %>%
  filter(k_mean > 0) %>%
  mutate(light_cat = factor(light_mean)) %>% # adding column for coloration purposes
  ggplot(aes(x = light_mean, y = k_mean, fill = light_cat)) +
  geom_point(shape = 21, size = 5, alpha = 0.75) +
  scale_fill_manual(values = cal_palette("canary", n = 28, type = "continuous")) + # custom colors
  labs(x = "Light Availability (ppfd?)",
       y = "Carrying Capacity (K)") +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")

fig_light2

#### Divergences ####

# fig_div <- ggplot(data_out_divs_df, aes(x = bin_cat, fill = bin_cat)) + # base plot
#   geom_histogram(stat = "count") + # divergences histogram
#   scale_fill_manual(values = cal_palette("fire", n = 7, type = "continuous")) + # custom colors
#   labs(x = "Number of Divergences",
#        y = "Site Count") +
#   scale_y_continuous(breaks=seq(0, 10, 1)) + # fix y axis labels
#   theme_classic() + # remove grid
#   theme(legend.position = "none")
# 
# fig_div
# DONT USE THIS FIGURE YET - for whatever reason, the 
# number of divergences on the shinystan app does not
# match the output of the above function from the mc-stan
# help document....

# For the meeting on 9/27/2021, I will be using the output of shinystan since,
# it's the larger/more conservative

data_divs_join <- left_join(data_years, data_out_diff_divs, by = "site_name")

# Adding category to make a nicer looking plots
data_divs_join <- data_divs_join %>%
  mutate(bin_cat = factor(case_when(div_shinyStan == 0 ~ "0",
                                    div_shinyStan > 0 & div_shinyStan <= 50 ~ "0-50",
                                    div_shinyStan > 50 & div_shinyStan <= 100 ~ "50-100",
                                    div_shinyStan > 100 & div_shinyStan <= 1000 ~ "100-1000",
                                    div_shinyStan > 1000 & div_shinyStan <= 2500 ~ "1000-2500",
                                    div_shinyStan > 2500 ~ "2500+",
                                    TRUE ~ "NA"),
                          levels = c("0", "0-50", "50-100", "100-1000", "1000-2500", "2500+")))

fig_yrs_divs <- ggplot(data_divs_join,
                       aes(x = n, y = div_shinyStan, fill = bin_cat, label = site_name)) +
  geom_point(shape = 21, size = 5, alpha = 0.75) +
  scale_fill_manual(values = cal_palette("fire", n = 6, type = "continuous")) + # custom colors
  labs(x = "Years of Available Data",
       y = "Divergent Transitions") +
  geom_text_repel(data = subset(data_divs_join, n > 5 & div_shinyStan > 2000), size = 4) +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")

fig_yrs_divs

# Also creating figure of hydrology vs. divergences

# Calculate coefficient of variation (sd/mean) for each site
test_site <- data_in$nwis_01608500
test_sd <- sd(test_site$Q)
test_mean <- mean(test_site$Q)
test_cv <- test_sd/test_mean

# Now, to create a function...
calc_coeff_var <- function(df){
  sd <- sd(df$Q)
  mean <- mean(df$Q)
  sd/mean
}

# And now map this to the entire output list.
data_cvs <- map(data_in, calc_coeff_var)
# yipee

# And create a dataframe
data_cvs_divs <- map_df(data_cvs, ~as.data.frame(.x), .id="site_name") %>%
  rename(coeff_var = `.x`) %>%
  left_join(data_out_diff_divs, by = "site_name")

# Adding same categories as above to make a nicer looking plots
data_cvs_divs <- data_cvs_divs %>%
  mutate(bin_cat = factor(case_when(div_shinyStan == 0 ~ "0",
                                    div_shinyStan > 0 & div_shinyStan <= 50 ~ "0-50",
                                    div_shinyStan > 50 & div_shinyStan <= 100 ~ "50-100",
                                    div_shinyStan > 100 & div_shinyStan <= 1000 ~ "100-1000",
                                    div_shinyStan > 1000 & div_shinyStan <= 2500 ~ "1000-2500",
                                    div_shinyStan > 2500 ~ "2500+",
                                    TRUE ~ "NA"),
                          levels = c("0", "0-50", "50-100", "100-1000", "1000-2500", "2500+")))

fig_cv_divs <- ggplot(data_cvs_divs,
                       aes(x = coeff_var, y = div_shinyStan, fill = bin_cat, label = site_name)) +
  geom_point(shape = 21, size = 5, alpha = 0.75) +
  scale_fill_manual(values = cal_palette("fire", n = 6, type = "continuous")) + # custom colors
  labs(x = "Coefficient of Variation in Discharge",
       y = "Divergent Transitions") +
  geom_text_repel(data = subset(data_cvs_divs, coeff_var > 3 | div_shinyStan > 3000), size = 4) +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")

fig_cv_divs

# creating composite figure for export
fig_div_full <- fig_yrs_divs + fig_cv_divs +
  plot_annotation(tag_levels = 'A')

fig_div_full

# ggsave(fig_div_full,
#        filename = "figures/teton_34sites/divergences_paneled.jpg",
#        width = 12,
#        height = 5)

#### Check model diagnostics ####

# Mimic the s vs. c plot created above, but this time with confidence intervals added.
# Careful, I am not filtering out negative r and k values here for now (just for a rough
# estimate of intervals about the mean).
data_summary_wide <- data_summary_siteinfo %>%
  filter(parameter == "s" | parameter == "c") %>%
  select(site_name, parameter, mean, `2.5%`, `97.5%`) %>%
  pivot_wider(names_from = parameter, values_from = c(mean, `2.5%`, `97.5%`))

fig_sc_ci.1 <- ggplot(data_summary_wide, aes(x = mean_c, y = mean_s, 
                                   fill = site_name)) + # , label = site_name
  geom_point(shape = 21, size = 4, alpha = 0.75) +
  labs(x = "Critical Discharge (c)",
       y = "Sensitivity of Persistence Curve (s)") +
  xlim(0, 3) +
  ylim(0, 600) +
  geom_errorbar(aes(ymin = `2.5%_s`,ymax = `97.5%_s`)) + 
  #geom_errorbarh(aes(xmin = `2.5%_c`,xmax = `97.5%_c`)) +
  #geom_text_repel(size=3) +
  scale_fill_manual(values = cal_palette("creek", n = 34, type = "continuous")) + # custom colors
  theme_bw() +
  theme(legend.position = "none")

fig_sc_ci.1

fig_sc_ci.2 <- ggplot(data_summary_wide, aes(x = mean_c, y = mean_s, 
                                             fill = site_name)) + # , label = site_name
  geom_point(shape = 21, size = 4, alpha = 0.75) +
  labs(x = "Critical Discharge (c)",
       y = "Sensitivity of Persistence Curve (s)") +
  xlim(0, 3) +
  ylim(0, 600) +
  #geom_errorbar(aes(ymin = `2.5%_s`,ymax = `97.5%_s`)) + 
  geom_errorbarh(aes(xmin = `2.5%_c`,xmax = `97.5%_c`)) +
  #geom_text_repel(size=3) +
  scale_fill_manual(values = cal_palette("creek", n = 34, type = "continuous")) + # custom colors
  theme_bw() +
  theme(legend.position = "none")

fig_sc_ci.2

fig_sc_ci_full <- fig_sc_ci.1 + fig_sc_ci.2

fig_sc_ci_full

# Export figure, but don't include in RMarkdown for clarity's sake.
# ggsave(fig_sc_ci_full,
#        filename = "figures/teton_34sites/s_vs_c_withcis.jpg",
#        width = 12,
#        height = 5)

# End of script.
