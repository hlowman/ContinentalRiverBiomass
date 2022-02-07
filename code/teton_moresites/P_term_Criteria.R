## P term inclusion criteria for 207 rivers dataset
## Created: February 7, 2022
## Heili Lowman

# The following code will filter the Teton output for the 207 site dataset,
# and based on a set of decided upon criteria, sites will be filtered out
# leaving behind the ones that should be re-run on Teton using a model structure
# that doesn't include a P term.

# Load packages
lapply(c("calecopal", "cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", "here",
         "ggrepel", "patchwork", "grid","gridExtra",
         "rlist", "pipeR"), require, character.only=T)

#### Data Import ####

# Load dataset loaded into Teton.
data_in <- readRDS("data_working/df_207sites.rds")

# Load dataset with site information.
data_info <- readRDS("data_working/NWIS_207sitesinfo_subset.rds")

# Load dataset with parameters' data (all iterations).
data_params <- readRDS("data_working/teton_207rivers_model_parameters_all_iterations_020622.rds")

# Load dataset with model divergences for each site.
data_divs <- readRDS("data_working/teton_207rivers_model_divergences_020622.rds")

# Load dataset with model diagnostics for each site.
data_diags <- readRDS("data_working/teton_207rivers_model_diagnostics_020622.rds")

#### Filtering ####

# Criteria:

# (1) If the Rhat values for c and s parameters are above 1.05, re-run.
c1 <- data_diags %>%
  filter(parameter == "c") %>%
  filter(Rhat > 1.05) %>% # 42 sites
  select(site_name)

s1 <- data_diags %>%
  filter(parameter == "s") %>%
  filter(Rhat > 1.05) %>% # 21 sites
  select(site_name)

# (2) If the n_eff value for c and s parameters is below 10% of
# the total numer of iterations (7500), re-run.
c2 <- data_diags %>%
  filter(parameter == "c") %>%
  filter(n_eff < 750) %>% # 103 sites
  select(site_name)

s2 <- data_diags %>%
  filter(parameter == "s") %>%
  filter(n_eff < 750) %>% # 49 sites
  select(site_name)

# (3) If the overall model fit had more than 50 divergences, re-run.

d3 <- data_divs %>%
  filter(divergences > 50) %>% # 81 sites
  select(site_name)

# Now, join together the lists of site names above to determine the sites
# that will be re-run.

join1 <- full_join(c1, s1) # 42 sites
join2 <- full_join(join1, c2) # 103 sites
join3 <- full_join(join2, s2) # 109 sites
joined <- full_join(join3, d3) # 132 sites

# And finally, to filter the data list generated earlier for the sites in
# "joined".

# Renaming the sites to be re-run just for clarity.
rerun <- joined

# Create vector of sites that will keep the P term in.
all_sites <- data_divs %>%
  select(site_name)

keep <- anti_join(all_sites, rerun)

# Remove these sites from the original dataset.
data_in_rerun <- list.remove(data_in, c("nwis_01124000", "nwis_01400500", 
"nwis_01408500", "nwis_01463500", "nwis_01480617", "nwis_01567000", 
"nwis_01608500", "nwis_01632900", "nwis_01646305", "nwis_01648010",
"nwis_01650800", "nwis_01654500", "nwis_01656903", "nwis_01673000",
"nwis_02110500", "nwis_02135000", "nwis_02160105", "nwis_02169000",
"nwis_02203950", "nwis_02204037", "nwis_02207135", "nwis_02207160",
"nwis_02208450", "nwis_02208493", "nwis_02266200", "nwis_02336120",
"nwis_02336340", "nwis_02336360", "nwis_02336526", "nwis_02344630",
"nwis_03025500", "nwis_03036000", "nwis_03061000", "nwis_03067510",
"nwis_03098600", "nwis_03106000", "nwis_03183500", "nwis_03292475",
"nwis_03292480", "nwis_03353200", "nwis_03353420", "nwis_03408500",
"nwis_03538830", "nwis_04108660", "nwis_04121944", "nwis_04124200",
"nwis_04136500", "nwis_04137005", "nwis_04137500", "nwis_041482663",
"nwis_04166500", "nwis_04167150", "nwis_04176500", "nwis_04199500",
"nwis_04200500", "nwis_04208000", "nwis_05054000", "nwis_05057000",
"nwis_05435943", "nwis_07048600", "nwis_07075250", "nwis_07143672",
"nwis_07144100", "nwis_07191222", "nwis_08057410", "nwis_08070200",
"nwis_08171290", "nwis_08180700", "nwis_11044000", "nwis_11462500",
"nwis_11463000", "nwis_11463980", "nwis_14301500", "nwis_385446094430700",
"nwis_385520094420000")) # 132 elements in the new list

# Export for run on Teton.
saveRDS(data_in_rerun, file = "data_working/df_207sites_rerun132.rds")

# End of script.
